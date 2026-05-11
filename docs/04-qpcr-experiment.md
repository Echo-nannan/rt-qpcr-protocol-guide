# 04. qPCR Experiment

## Goal

Generate reliable Ct values with a stable reaction system and plate layout, using NTC, NRT, technical replicates, and melt curves for QC.

## Example qPCR Reaction Setup

### 20 uL Reaction

| Component | Volume |
|---|---:|
| 2x SYBR Green qPCR Mix | 10 uL |
| Forward primer, 10 uM | 0.4 uL |
| Reverse primer, 10 uM | 0.4 uL |
| cDNA template | 2 uL |
| RNase-free water | 7.2 uL |
| Total | 20 uL |

### 10 uL Reaction

| Component | Volume |
|---|---:|
| 2x SYBR Green qPCR Mix | 5 uL |
| Forward primer, 10 uM | 0.4 uL |
| Reverse primer, 10 uM | 0.4 uL |
| cDNA template | 1 uL |
| RNase-free water | 3.2 uL |
| Total | 10 uL |

## Example Cycling Program

| Step | Temperature | Time | Cycles |
|---|---:|---:|---:|
| Initial denaturation | 95 degC | 30 s | 1 |
| Denaturation | 95 degC | 10 s | 40 |
| Annealing/extension | 60 degC | 30 s | 40 |
| Melt curve | 60-95 degC | ramp | 1 |

## Plate Layout Rules

- Use three technical replicates for each sample-gene combination.
- Include at least one to three NTC wells for each primer pair.
- Include NRT controls to check genomic DNA contamination.
- Reduce edge effects; add water or buffer to edge wells when needed.
- Group samples or genes consistently to reduce pipetting mistakes.

## Post-Run Checklist

- [ ] Amplification curves show standard sigmoidal behavior.
- [ ] Melt curves show a single peak.
- [ ] NTC has no amplification or Ct >= 38.
- [ ] NRT has no amplification.
- [ ] Technical replicate SD is < 0.5 Ct.
- [ ] Reference gene Ct is stable across samples.
