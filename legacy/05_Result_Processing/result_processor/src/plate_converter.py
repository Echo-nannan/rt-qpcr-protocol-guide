"""
PCR数据转换为384孔板布局模块
"""
import re
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd


class PCRPlateConverter:
    """将qPCR结果文件转换为384孔板布局的Excel格式"""
    
    ROWS = list('ABCDEFGHIJKLMNOP')  # 16行
    COLS = list(range(1, 25))  # 24列
    
    # 样本配对: A+B, C+D, E+F, G+H, I+J, K+L, M+N, O+P
    ROW_PAIRS = [('A', 'B'), ('C', 'D'), ('E', 'F'), ('G', 'H'),
                 ('I', 'J'), ('K', 'L'), ('M', 'N'), ('O', 'P')]
    
    def __init__(self):
        self.data: List[Dict] = []
    
    def parse_lightcycler(self, filepath: str) -> List[Dict]:
        """解析Roche LightCycler导出的文件"""
        data = []
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        for line in lines[2:]:  # 跳过实验信息和表头
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                pos = parts[2]
                match = re.match(r'([A-P])(\d+)', pos)
                if match:
                    data.append({
                        'Include': parts[0],
                        'Color': parts[1],
                        'Pos': pos,
                        'Name': parts[3],
                        'Cp': parts[4].strip() if parts[4] else '',
                        'Row': match.group(1),
                        'Col': int(match.group(2))
                    })
        self.data = data
        return data
    
    def to_plate_layout(self, output_path: str) -> None:
        """导出为384孔板布局的Excel文件（原始格式）"""
        if not self.data:
            raise ValueError("没有数据，请先调用parse方法")
        
        data_dict = {(d['Row'], d['Col']): d for d in self.data}
        excel_data = []
        
        for row in self.ROWS:
            row_data = []
            for col in self.COLS:
                key = (row, col)
                if key in data_dict:
                    d = data_dict[key]
                    row_data.extend([
                        d['Include'], d['Color'], d['Pos'], 
                        d['Name'], d['Cp'], '', '0', ''
                    ])
                else:
                    row_data.extend([''] * 8)
            excel_data.append(row_data)
        
        columns = []
        for col in self.COLS:
            columns.extend([
                f'Include_{col}', f'Color_{col}', f'Pos_{col}', f'Name_{col}',
                f'Cp_{col}', f'Conc_{col}', f'Std_{col}', f'Sep_{col}'
            ])
        
        df = pd.DataFrame(excel_data, columns=columns, index=self.ROWS)
        
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Plate_Layout', index=True)
        
        print(f"✅ 已保存: {output_path}")
    
    def to_paired_layout(self, output_path: str) -> None:
        """
        导出为配对布局的Excel文件
        每两行(A+B, C+D等)为一个样本，按列1-24竖着排列
        便于直接复制粘贴
        """
        if not self.data:
            raise ValueError("没有数据，请先调用parse方法")
        
        data_dict = {(d['Row'], d['Col']): d for d in self.data}
        
        # 每个样本对生成一组列
        all_columns = []
        
        for pair_idx, (row1, row2) in enumerate(self.ROW_PAIRS):
            # 每个配对样本的数据：24行（对应1-24列），每行包含两个孔的数据
            pair_data = []
            for col in self.COLS:
                # 第一行数据 (如A1)
                key1 = (row1, col)
                d1 = data_dict.get(key1, {})
                # 第二行数据 (如B1)
                key2 = (row2, col)
                d2 = data_dict.get(key2, {})
                
                pair_data.append({
                    'col': col,
                    # 第一行
                    f'{row1}_Include': d1.get('Include', ''),
                    f'{row1}_Color': d1.get('Color', ''),
                    f'{row1}_Pos': d1.get('Pos', ''),
                    f'{row1}_Name': d1.get('Name', ''),
                    f'{row1}_Cp': d1.get('Cp', ''),
                    f'{row1}_Std': '0' if d1 else '',
                    # 第二行
                    f'{row2}_Include': d2.get('Include', ''),
                    f'{row2}_Color': d2.get('Color', ''),
                    f'{row2}_Pos': d2.get('Pos', ''),
                    f'{row2}_Name': d2.get('Name', ''),
                    f'{row2}_Cp': d2.get('Cp', ''),
                    f'{row2}_Std': '0' if d2 else '',
                })
            
            all_columns.append((f'{row1}+{row2}', pair_data))
        
        # 创建Excel，每个样本对占一组列
        excel_rows = []
        for col_idx in range(24):  # 24行对应1-24列
            row_data = []
            for pair_name, pair_data in all_columns:
                d = pair_data[col_idx]
                row1, row2 = pair_name.split('+')
                row_data.extend([
                    d[f'{row1}_Include'], d[f'{row1}_Color'], d[f'{row1}_Pos'],
                    d[f'{row1}_Name'], d[f'{row1}_Cp'], '', d[f'{row1}_Std'], '',
                    d[f'{row2}_Include'], d[f'{row2}_Color'], d[f'{row2}_Pos'],
                    d[f'{row2}_Name'], d[f'{row2}_Cp'], '', d[f'{row2}_Std'], ''
                ])
            excel_rows.append(row_data)
        
        # 创建列名
        columns = []
        for pair_name, _ in all_columns:
            row1, row2 = pair_name.split('+')
            for row in [row1, row2]:
                columns.extend([
                    f'{row}_Inc', f'{row}_Color', f'{row}_Pos', f'{row}_Name',
                    f'{row}_Cp', f'{row}_Conc', f'{row}_Std', ''
                ])
        
        df = pd.DataFrame(excel_rows, columns=columns)
        df.index = self.COLS  # 行索引为1-24
        
        with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Paired_Layout', index=True)
        
        print(f"✅ 已保存: {output_path}")
    
    def convert(self, input_path: str, output_path: Optional[str] = None, 
                paired: bool = True) -> str:
        """
        一键转换文件
        
        Args:
            input_path: 输入文件路径
            output_path: 输出文件路径
            paired: True使用配对布局(A+B为一组)，False使用原始布局
        """
        input_path = Path(input_path)
        if output_path is None:
            output_path = input_path.parent / f"{input_path.stem}_plate.xlsx"
        
        self.parse_lightcycler(str(input_path))
        
        if paired:
            self.to_paired_layout(str(output_path))
        else:
            self.to_plate_layout(str(output_path))
        
        return str(output_path)
