# 06. Automation Pipeline

## 目标

把人工整理 Ct 表、质控、Delta Delta Ct 计算、统计检验和图表导出整合成可重复运行的分析流程。

## 推荐流水线

```mermaid
flowchart LR
    A["输入文件"] --> B["Schema 校验"]
    B --> C["Ct 质控"]
    C --> D["Delta Ct"]
    D --> E["Delta Delta Ct"]
    E --> F["统计与 FDR"]
    F --> G["图表与报告"]
    G --> H["归档 provenance"]
```

## 输入契约

### `sample_info.csv`

| 列名 | 类型 | 说明 |
|---|---|---|
| sample_id | string | 样本唯一 ID |
| group | string | 分组，如 Control/Treatment |
| tissue | string | 可选，组织或细胞类型 |
| batch | string | 可选，实验批次 |

### `ct_values.csv`

| 列名 | 类型 | 说明 |
|---|---|---|
| sample_id | string | 对应样本 ID |
| gene | string | 目标基因或内参基因 |
| ct | number | Ct 值 |
| replicate | integer | 技术重复编号 |
| well | string | 可选，孔位 |

## 输出目录建议

```text
results/
└── run_YYYYMMDD_HHMMSS/
    ├── tables/
    │   ├── detail.csv
    │   └── summary.csv
    ├── figures/
    │   ├── fold_change.png
    │   └── volcano.png
    ├── report.md
    └── provenance.json
```

## 自动化要点

- 明确记录 reference gene、control group、统计方法和离群值方法。
- 所有输入输出保存相对路径，方便迁移。
- 每次运行生成唯一 run id。
- 保留原始 Ct 表，不在原文件上修改。
- 图表和表格由同一份清洗后数据生成，避免结果不一致。

