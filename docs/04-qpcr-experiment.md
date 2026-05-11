# 04. qPCR Experiment

## 目标

通过稳定的反应体系和板布局获得可靠 Ct 值，并用 NTC、NRT、技术重复和熔解曲线进行质量控制。

## qPCR 反应体系示例

### 20 uL 体系

| 组分 | 体积 |
|---|---:|
| 2x SYBR Green qPCR Mix | 10 uL |
| Forward primer, 10 uM | 0.4 uL |
| Reverse primer, 10 uM | 0.4 uL |
| cDNA 模板 | 2 uL |
| RNase-free water | 7.2 uL |
| 总计 | 20 uL |

### 10 uL 体系

| 组分 | 体积 |
|---|---:|
| 2x SYBR Green qPCR Mix | 5 uL |
| Forward primer, 10 uM | 0.4 uL |
| Reverse primer, 10 uM | 0.4 uL |
| cDNA 模板 | 1 uL |
| RNase-free water | 3.2 uL |
| 总计 | 10 uL |

## 程序示例

| 阶段 | 温度 | 时间 | 循环 |
|---|---:|---:|---:|
| 预变性 | 95 degC | 30 s | 1 |
| 变性 | 95 degC | 10 s | 40 |
| 退火/延伸 | 60 degC | 30 s | 40 |
| 熔解曲线 | 60-95 degC | 递增 | 1 |

## 板布局原则

- 每个样本-基因组合设置 3 个技术重复。
- NTC 每个引物对至少设置 1-3 孔。
- NRT 用于检查 gDNA 污染。
- 尽量减少边缘效应，必要时边缘孔加水或 buffer。
- 同一基因或同一样本尽量成组排布，降低移液错误。

## 上机后检查

- [ ] 扩增曲线为标准 S 型
- [ ] 熔解曲线为单一峰
- [ ] NTC 无扩增或 Ct >= 38
- [ ] NRT 无扩增
- [ ] 技术重复 SD < 0.5 Ct
- [ ] 内参基因 Ct 在各样本间稳定

