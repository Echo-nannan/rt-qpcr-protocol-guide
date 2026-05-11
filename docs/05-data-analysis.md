# 05. Data Analysis

## 输入数据

建议把仪器导出的 Ct 值整理成长表格式：

| sample_id | gene | ct | replicate | well |
|---|---|---:|---:|---|
| S1 | GAPDH | 18.5 | 1 | A1 |
| S1 | GeneA | 24.2 | 1 | B1 |

样本分组表：

| sample_id | group | tissue | batch |
|---|---|---|---|
| S1 | Control | Brain | Batch1 |
| S2 | Treatment | Brain | Batch1 |

## Delta Delta Ct 计算

```text
Delta Ct = Ct_target - Ct_reference

Delta Delta Ct = Delta Ct_treatment - mean(Delta Ct_control)

Fold Change = 2^(-Delta Delta Ct)

log2FC = -Delta Delta Ct
```

## 多内参归一化

多内参时推荐使用内参 Ct 的几何均值思想进行归一化，避免单一内参波动造成偏差。

```text
reference Ct = geometric mean of reference gene Cts
Delta Ct = target Ct - reference Ct
```

## 统计建议

| 场景 | 推荐方法 |
|---|---|
| 两组比较，方差近似 | Student t-test |
| 两组比较，方差不齐 | Welch t-test |
| 两组非正态 | Mann-Whitney U |
| 多组比较 | one-way ANOVA |
| 多基因多重检验 | Benjamini-Hochberg FDR |

统计检验优先在 Delta Ct 上进行，而不是直接在 fold change 上做检验。

## 显著性标记

| 标记 | 阈值 |
|---|---:|
| ns | p 或 FDR >= 0.05 |
| * | < 0.05 |
| ** | < 0.01 |
| *** | < 0.001 |

## 输出建议

- `summary.csv`: 每个基因的 fold change、log2FC、p value、FDR
- `detail.csv`: 每个样本的 Ct、Delta Ct、Delta Delta Ct
- `figures/`: 柱状图、散点图、火山图
- `provenance.json`: 参数、输入文件 hash、软件版本和运行时间

