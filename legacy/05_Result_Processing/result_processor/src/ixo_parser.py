"""
Roche LightCycler 480 (.ixo) 文件解析器
========================================

.ixo 是 Roche LightCycler 480 软件 (LCS480) 的原生实验文件，本质是一种
名为 IXOS 的对象序列化文本格式（伪 XML）。它包含了实验全部信息：
metadata、原始荧光曲线、分析结果（CrossingPoint / Cp 值）等。

本模块只关注 384 个孔的 Cp（CT）值与孔位映射，输出格式与
``PCRDataParser.parse_file()`` 完全一致，可直接复用现有 ΔΔCt 流水线。

文件结构关键片段
----------------
.. code-block:: xml

    <list name="AnaSamples" count="384">
      <obj name="item" class="QuantSampleB" version="2">
        <prop name="IsIncluded">1</prop>
        <prop name="Pos">0</prop>            <!-- 0..383 行优先 -->
        <prop name="CrossingPoint">15.6647637787881</prop>
        ...
      </obj>
      ...
    </list>

孔位映射规则::

    Pos  = row_index * 24 + col_index   (0-indexed, row-major)
    Pos  = 0   -> A1
    Pos  = 23  -> A24
    Pos  = 24  -> B1
    Pos  = 383 -> P24
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional


ROWS = list("ABCDEFGHIJKLMNOP")  # 16 行 A-P


def pos_to_well(pos: int) -> tuple[str, int]:
    """0-indexed Pos → (row letter, 1-indexed col)。"""
    if pos < 0 or pos > 383:
        raise ValueError(f"Pos 越界: {pos}（应在 0..383）")
    row_idx, col_offset = divmod(pos, 24)
    return ROWS[row_idx], col_offset + 1


def is_ixo_file(filepath: str | Path) -> bool:
    """通过签名判断是否为 .ixo 文件（容忍后缀大小写或被改名）。"""
    p = Path(filepath)
    if p.suffix.lower() == ".ixo":
        return True
    try:
        with open(p, "rb") as f:
            head = f.read(64)
        return b'signature="IXOS"' in head or head.lstrip().startswith(b"<objectstream")
    except OSError:
        return False


class IxoParser:
    """解析 .ixo 文件并提取每孔 CrossingPoint。

    用法::

        parser = IxoParser()
        parser.parse_file("xxx.ixo")
        # parser.data 与 PCRDataParser.data 字段完全一致

    Attributes
    ----------
    data : List[Dict]
        每个 dict 含 ``Pos``(如 ``"A1"``)、``Name``(默认占位)、``Cp``(字符串)、
        ``Row``、``Col``。Cp 为空字符串表示该孔未检出（CrossingPoint 缺失）。
    experiment_name : str
        从 ``<prop name="name">`` 提取到的实验名（用于显示，可选）。
    """

    # 匹配单个 QuantSampleB 块。class 必须是 QuantSampleB 才认（避免误匹配
    # MeltingSample 等其他节点）。
    _BLOCK_RE = re.compile(
        r'<obj\s+name="item"\s+class="QuantSampleB"[^>]*>(.*?)</obj>',
        re.DOTALL,
    )
    _POS_RE = re.compile(r'<prop\s+name="Pos">\s*(-?\d+)\s*</prop>')
    _CP_RE = re.compile(
        r'<prop\s+name="CrossingPoint">\s*([0-9.+\-eE]+)\s*</prop>'
    )
    _CP_STATUS_RE = re.compile(
        r'<prop\s+name="CrossingPointStatus">\s*(-?\d+)\s*</prop>'
    )
    _CALL_RE = re.compile(
        r'<prop\s+name="Call">\s*(\d+)\s*</prop>'
    )
    _IS_INCLUDED_RE = re.compile(
        r'<prop\s+name="IsIncluded">\s*(\d+)\s*</prop>'
    )
    _EXP_NAME_RE = re.compile(
        r'<obj\s+name="root"[^>]*>\s*<prop\s+name="name">\s*([^<]*?)\s*</prop>',
        re.DOTALL,
    )

    def __init__(self) -> None:
        self.data: List[Dict] = []
        self.data_dict: Dict[tuple, Dict] = {}
        self.experiment_name: str = ""
        self._raw_text: str = ""

    # ------------------------------------------------------------------ public
    def parse_file(self, filepath: str | Path) -> List[Dict]:
        """解析 .ixo 文件，返回 384 条记录。"""
        p = Path(filepath)
        # IXOS 大体是文本（latin-1 即可全字节穿透），里面偶尔混入二进制 blob，
        # 但 Pos / CrossingPoint 始终是纯文本，所以 latin-1 解码最稳妥。
        with open(p, "rb") as f:
            raw = f.read()
        text = raw.decode("latin-1", errors="replace")
        self._raw_text = text

        m = self._EXP_NAME_RE.search(text)
        if m:
            self.experiment_name = m.group(1).strip()

        records: Dict[int, Dict] = {}
        for block_match in self._BLOCK_RE.finditer(text):
            block = block_match.group(1)
            pos_m = self._POS_RE.search(block)
            if not pos_m:
                continue
            pos = int(pos_m.group(1))
            if pos < 0 or pos > 383:
                continue
            # 以 IsIncluded 过滤（Roche 标记为 0 的孔表示用户取消），但默认仍写入
            # 数据（Cp 留空即可），保持 384 完整布局。
            included_m = self._IS_INCLUDED_RE.search(block)
            included = (
                included_m is None or included_m.group(1) == "1"
            )

            cp_m = self._CP_RE.search(block)
            call_m = self._CALL_RE.search(block)
            cp_str = ""
            if cp_m and included:
                # Roche Call 字段语义：
                #   Call=2  完全有效 Cp                     → 保留
                #   Call=1  Detector Call uncertain (有值)  → 保留（与 Roche TXT 导出一致）
                #   Call=0  未检出（Cp 字段填占位 0）        → 留空
                # 仅以 Call=0 为"未检出"判据；其余情况一律取 CrossingPoint。
                call_value = call_m.group(1) if call_m else "2"
                if call_value != "0":
                    try:
                        cp_val = float(cp_m.group(1))
                        # 真实 Cp 不可能为 0；如果 Roche 给 0 即视为未检出。
                        if cp_val > 0:
                            cp_str = f"{cp_val:.2f}"
                    except ValueError:
                        cp_str = ""

            row_letter, col = pos_to_well(pos)
            records[pos] = {
                "Pos": f"{row_letter}{col}",
                "Name": f"Sample {pos + 1}",
                "Cp": cp_str,
                "Row": row_letter,
                "Col": col,
            }

        # 对未出现的孔位补空白记录，保证下游遍历 384 孔时不会缺位。
        data: List[Dict] = []
        for pos in range(384):
            if pos in records:
                data.append(records[pos])
            else:
                row_letter, col = pos_to_well(pos)
                data.append(
                    {
                        "Pos": f"{row_letter}{col}",
                        "Name": f"Sample {pos + 1}",
                        "Cp": "",
                        "Row": row_letter,
                        "Col": col,
                    }
                )

        self.data = data
        self.data_dict = {(d["Row"], d["Col"]): d for d in data}
        return data

    def to_lightcycler_text(self, output_path: str | Path) -> Path:
        """将 .ixo 数据导出为 Roche 标准文本格式（与 .txt 导出兼容）。

        便于下游沿用现有 ``PCRDataParser`` 的解析逻辑做交叉验证。
        """
        if not self.data:
            raise ValueError("数据为空，请先调用 parse_file()")
        out = Path(output_path)
        with open(out, "w", encoding="utf-8") as f:
            f.write(
                f"Experiment: {self.experiment_name or 'Imported from IXO'}\n"
            )
            f.write("Include\tColor\tPos\tName\tCp\tConcentration\tStandard\tStatus\n")
            for d in self.data:
                f.write(
                    "True\t255\t"
                    f"{d['Pos']}\t{d['Name']}\t{d['Cp']}\t\t0\t\n"
                )
        return out


__all__ = ["IxoParser", "is_ixo_file", "pos_to_well", "ROWS"]
