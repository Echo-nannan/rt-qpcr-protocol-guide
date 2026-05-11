# 03. Reverse Transcription

## 目标

将等量 RNA 逆转录为 cDNA，保证不同样本在 qPCR 阶段具有可比性。

## 计算原则

RNA 投入量建议在同一批实验中保持一致。例如统一投入 500 ng 或 1000 ng。

```text
RNA 体积 uL = 目标投入量 ng / RNA 浓度 ng/uL
```

如果计算得到 RNA 体积过大，应考虑浓缩 RNA 或降低所有样本的统一投入量。

## 20 uL RT 体系示例

| 组分 | 体积 | 说明 |
|---|---:|---|
| RNA 模板 | X uL | 含等量总 RNA |
| 4x gDNA wiper Mix | 4 uL | 去除 gDNA |
| RNase-free water | 补至 16 uL | 第一阶段体系 |
| 5x qRT SuperMix | 4 uL | 逆转录体系 |
| 总体积 | 20 uL |  |

## 程序示例

| 阶段 | 温度 | 时间 |
|---|---:|---:|
| gDNA 去除 | 42 degC | 2 min |
| 引物退火 | 25 degC | 5 min |
| 逆转录 | 50 degC | 15 min |
| 酶灭活 | 85 degC | 5 sec |
| 保存 | 4 degC | hold |

## cDNA 使用

- 短期可 4 degC 保存，长期建议 -20 degC。
- qPCR 前常用 5x 或 10x 稀释。
- 同一实验所有样本使用相同稀释倍数。
- 避免反复冻融，必要时分装保存。

