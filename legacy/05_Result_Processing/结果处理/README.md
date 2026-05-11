# qPCR结果处理器 Pro

🧬 一款功能强大的qPCR数据处理工具，支持格式转换、ΔΔCt计算和自定义导出。

## 功能特性

### 📋 格式转换
- 解析Roche LightCycler等仪器导出的原始数据文件
- 支持384孔板布局（24列 × 16行）
- 可自定义样本名、基因名、分组名
- 导出原始数据、平均CT值、长格式数据（适合ΔΔCt计算）

### 🧮 ΔΔCt计算
- 完整的ΔΔCt相对定量计算
- 支持按组计算内参均值
- 自动检测基因和样本列表
- 输出详细计算过程和汇总表
- **并排格式导出**：内参基因和目的基因对比显示

### 📝 自定义输出文件名
支持模板变量：
- `{date}` - 当前日期（格式：YYYYMMDD）
- `{sample}` - 样本信息（可自定义）
- `{gene}` - 基因名

**示例**：`{date}_{sample}_Analysis and processing` → `20250124_DemoSample_Analysis and processing.csv`

## 安装

### Python环境运行
```bash
# 安装依赖
pip install -r requirements.txt

# 运行
python scripts/start_gui.py

# 或者
run_gui.bat
```

## 依赖
- Python 3.8+
- pandas >= 1.5.0
- openpyxl >= 3.0.0
- numpy >= 1.21.0
- xlsxwriter >= 3.0.0
- pyyaml >= 6.0.0
- matplotlib >= 3.7.0
- **scipy >= 1.10.0**（v2.4：t / Welch / Mann-Whitney 检验）
- **Pillow >= 9.0.0**（v2.4：openpyxl 嵌入 Fold-change PNG 图）
- tkinter（Python标准库）

## 使用方法

### 格式转换（支持 .txt 与 .ixo）
1. 点击"浏览文件"选择 qPCR 原始结果文件 — **既可选 Roche 文本导出 `.txt`，也可直接选 `.ixo` 原生实验文件**，自动识别。
2. 设置布局参数（每样本列数、每基因行数）
3. 可选：输入自定义样本名、基因名（逗号分隔）
4. 点击"预览"查看结果，或直接"导出 Excel"
5. 也可点击底部「🚀 一键 ΔΔCt 全流程」直接产出并排格式 ΔΔCt（弹窗选好内参基因 + 对照组即可）。

### 宽表 → ΔΔCt（v2.2 新增，最便捷的工作流）
当你已经有「Sample × Gene」宽格式 Excel/CSV（本工具导出的「原始数据」/「平均CT值」Sheet 都满足，或者从其他途径整理得到的也行）：
1. **Step 1**：选择宽格式 Excel / CSV（自动列出所有 Sheet 让你选）。
2. **Step 2**：自动识别样本/分组列与基因列，必要时手动选定样本列。
3. **Step 3**：下拉选**内参基因**（默认自动猜 GAPDH / ACTB / RPL13a 等常见名字）。
4. **Step 4**：在多选 checkbox 里勾选要分析的**目标基因**（默认全选除内参外）。
5. **Step 5**：下拉选**对照组**（自动从 Sample 列里抽出 `C / TD-L / TD-M / TD-H` 这种唯一标签）。
6. **Step 6**：（可选）调整文件名模板与样本信息。
7. **Step 7**：点击「🔢 计算 ΔΔCt」，然后选「📄 导出并排格式 CSV」或「💾 导出 Excel」。

### ΔΔCt 计算（长格式输入，原始流程保留）
1. **Step 1**：选择 Ct 数据文件（Excel/CSV）
   - 推荐使用格式转换生成的"长格式_ΔΔCt用"Sheet
2. **Step 2**：确认列名设置（通常自动检测）
3. **Step 3**：选择**主内参**（如 GAPDH）；可选 *多内参 (可选)* listbox 里 Ctrl+点击额外的内参基因（如同时选 ACTB / RPL13a）→ 几何平均合成综合内参。
4. **Step 4**：选择对照组
5. **Step 4·补**：选择**统计检验** (Student t / Welch t / Mann-Whitney U)、**图表格式** (PNG / PDF / SVG)、是否**嵌图到 Excel**；可点「💾 保存分析预设」一键存档便于下次复用。
6. **Step 5**：设置输出文件名模板
7. **Step 6**：点击"计算 ΔΔCt"，然后导出（柱状图自动加 ns/*/**/*** 显著性，Excel 自动嵌入图表 sheet）

### 输出格式

#### 并排格式CSV（推荐）
类似于：
```
Sample Name | Target Name | CT Value | 内参均值 | ... | Sample Name | Target Name | CT Value | ΔCt | ΔΔCt | 2^-ΔΔCt
    C6      |    GAPDH    |  16.69   |  15.64  | ... |     C6      |   Camk2a   |  15.61   | -0.1| -0.18|  1.13
```

#### Excel多Sheet格式
- **详细计算结果**：每个CT值一行，包含所有中间计算步骤
- **汇总表**：每个样本每个基因一行
- **并排格式_Analysis**：内参-目的基因对比
- **计算参数**：记录使用的设置

## 计算公式

```
① 组内参Ct均值 = 同一处理组内所有样本的内参基因Ct平均值
② 样本目的基因Ct均值 = 每个样本目的基因所有重复的Ct平均值
③ ΔCt = 样本目的基因Ct均值 - 组内参Ct均值
④ 对照组ΔCt均值 = 对照组所有样本的ΔCt平均值
⑤ ΔΔCt = 样本ΔCt - 对照组ΔCt均值
⑥ 2^(-ΔΔCt) = 相对表达量（Fold Change）
```

## 项目结构

```
结果处理/
├── presets/                   # 布局预设 JSON（用户可保存/加载）
├── src/
│   ├── complete_gui.py        # 主 GUI 程序（含 3 个标签页 + 一键全流程）
│   ├── ixo_parser.py          # Roche LightCycler 480 .ixo 原生文件解析器
│   └── plate_converter.py     # 板式转换模块
├── scripts/
│   └── start_gui.py           # 启动脚本
├── configs/
│   ├── example_layout.txt     # 布局配置示例
│   └── export_formats.yaml    # 导出格式配置
├── examples/
│   └── sample_ddct_input.csv  # 示例输入数据
├── requirements.txt
├── run_gui.bat
└── README.md
```

## 配置文件说明

### configs/example_layout.txt
```
sample_count=8
sample_0=Control1,1,3
sample_1=Control2,4,6
...
gene_count=4
gene_0=GAPDH,A,D
gene_1=GeneA,E,H
...
```

### configs/export_formats.yaml
可自定义导出格式、样本名、基因名和界面主题。

## 常见问题

**Q: 导出的CSV中文乱码？**
A: 程序使用UTF-8-BOM编码导出，Excel可正常识别。如仍有问题，请用Excel的"数据→从文本"导入，选择UTF-8编码。

**Q: 如何处理多个目的基因？**
A: 程序会自动识别所有非内参基因作为目的基因，分别计算并在并排格式中显示。

## 更新日志

### v2.5 (2026-04)
- **Pfaffl 法定量**：「ΔΔCt 计算」与「宽表 → ΔΔCt」两页都新增「定量方法」下拉（**`2^-ΔΔCt`** / **`Pfaffl (E^-ΔΔCt)`**）。选中 Pfaffl 后点 **「📊 设置基因效率…」** 弹出表格，为每个基因填扩增效率（接受 base `1.95` / 百分比 `95` / 增量 `0.85` 三种写法，会自动归一到 `1~2.5` 之间）；未填的基因回退到 base = 2.0（与 ΔΔCt 完全等价）。
- **效率表共享**：ΔΔCt 标签页 / 宽表标签页 / 一键全流程 三处共用同一个基因效率表，配置一次即可。
- **Excel 计算参数 sheet 增加** 「定量方法」 与 「基因效率 (Pfaffl)」 两行，方便溯源。
- 新 smoke 脚本 `scripts/_smoke_test_pfaffl.py`：默认 ΔΔCt 法 fold-change ≈ 2，Pfaffl(E=1.8) ≈ 1.8，Pfaffl(95%) ≈ 1.95 全部精确通过。

### v2.4 (2026-04)
- **多内参基因 + 几何平均**：「ΔΔCt 计算」与「宽表 → ΔΔCt」两页都新增「多内参 (可选)」listbox（``selectmode=MULTIPLE``）。多选时自动按 *几何平均* 合成一个综合内参 Ct（先在每组里按基因取均值再几何平均，避免不同基因量级互相挤压）；只勾「主内参」时与旧版完全等价。一键全流程的弹窗也同步加了多内参 + 统计方法 + 嵌图选项。
- **统计检验 + FDR + 显著性标星**：新增「统计检验」下拉（**Student t / Welch t / Mann-Whitney U**），ΔCt 对照组 vs 处理组检验 + Benjamini-Hochberg FDR；汇总表会自动多出 `n / fc_sem / pvalue / fdr / 显著性` 列。Fold-change 柱状图自动在每根柱顶标 `ns / * / ** / ***`。统计依赖 `scipy`（已加到 `requirements.txt`）。
- **图表导出 PNG / PDF / SVG**：新增「图表格式」下拉，可在三种格式间切换（PDF/SVG 直出矢量图，论文排版友好）。
- **导出 Excel 时嵌入柱状图**（可勾选）：所有「导出 Excel」入口（ΔΔCt、宽表、一键全流程）都会在 Excel 末尾追加一个 *Fold-change 图* sheet，直接嵌入 PNG（300 DPI）。需要 `Pillow + openpyxl`，未装时会在弹窗里给出原因，导出本身不会失败。
- **分析预设 (JSON)**：「ΔΔCt 计算」页新增「💾 保存 / 📂 加载分析预设」按钮，可把当前「主内参 / 多内参 / 对照组 / 统计方法 / 图表格式 / 嵌图开关 / 文件名模板」一键存到 `presets/`，下次回填。
- 老版 wide-tab 与 ΔΔCt 老接口完全兼容，旧的 `calculator.calculate(ref_gene=...)` 仍按单内参跑。

### v2.3 (2026-04)
- **布局预设**：「格式转换」页新增「布局预设 (JSON)」——一键保存/加载「每样本列数、每基因行数、每组样本数、样本/基因/处理组名称」，文件默认在 `presets/` 目录，方便同板型实验复用。
- **Fold-change 发表级柱状图**：在「ΔΔCt 计算」与「宽表 → ΔΔCt」页，计算完成后可选 **「(全部)」或单个目标基因**，导出 **PNG**（300 DPI，均值 ± SEM，虚线 y=1 参考）。依赖 `matplotlib`（已写入 `requirements.txt`）。
- **宽表质控**：「宽表 → ΔΔCt」Step 2 增加 **「缺失值统计」**，表格列出各基因列有效数值与空/非数个数。
- **CSV 自动嗅探**：宽表导入 CSV 时自动尝试 UTF-8 / GBK 编码与逗号、制表符、自动分隔符解析。

### v2.2 (2026-04)
- **新增 `.ixo` 直读支持**：Roche LightCycler 480 原生 IXOS 文件可直接喂给「格式转换」与「ΔΔCt 计算」标签页，无需先导出 .txt（解析与官方文本导出 100% 一致）。
- **新增「📊 宽表 → ΔΔCt」标签页**：直接吃 `Sample × Gene` 宽格式 Excel/CSV（如本工具导出的「原始数据」/「平均CT值」Sheet 或手工粘贴），通过下拉选**内参基因**、勾选**多个目标基因**、选择**对照组**，一键产出与「ΔΔCt 计算」同款的「内参/目的基因并排格式」CSV / 多 Sheet Excel。同分组下的多行会被自动赋予唯一样本 ID（如 `C_01..C_12`），保留原始分组语义。
- **新增「🚀 一键 ΔΔCt 全流程」按钮**（位于「格式转换」标签页底部）：选好 `.ixo` / `.txt` + 布局 + 自定义名称后，一个按钮 + 小弹窗即可输出含「长格式 / 详细 / 汇总 / 并排 / 参数」5 Sheet 的 Excel。
- 修复 `_load_config` 调用缺失的 `import yaml`，并增强 `.txt` 编码兼容（utf-8 / utf-8-sig / gbk / cp936 自动回退）。

### v2.1 (2025-01)
- 新增自定义输出文件名模板功能
- 新增并排格式CSV导出（内参-目的基因对比）
- 优化Excel导出，包含4个Sheet

### v2.0
- 完整GUI界面
- 支持格式转换和ΔΔCt计算
- 支持分组计算

## 许可证

MIT License
