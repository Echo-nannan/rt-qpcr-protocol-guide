# 05. Data Analysis

## Input Data

Organize instrument-exported Ct values in long-table format:

| sample_id | gene | ct | replicate | well |
|---|---|---:|---:|---|
| S1 | GAPDH | 18.5 | 1 | A1 |
| S1 | GeneA | 24.2 | 1 | B1 |

Sample metadata table:

| sample_id | group | tissue | batch |
|---|---|---|---|
| S1 | Control | Brain | Batch1 |
| S2 | Treatment | Brain | Batch1 |

## Delta Delta Ct Calculation

```text
Delta Ct = Ct_target - Ct_reference

Delta Delta Ct = Delta Ct_treatment - mean(Delta Ct_control)

Fold Change = 2^(-Delta Delta Ct)

log2FC = -Delta Delta Ct
```

## Multiple Reference Genes

When multiple reference genes are used, normalize with the geometric-mean concept for reference Ct values. This reduces bias from instability in a single reference gene.

```text
reference Ct = geometric mean of reference gene Cts
Delta Ct = target Ct - reference Ct
```

## Statistical Recommendations

| Scenario | Recommended method |
|---|---|
| Two groups with similar variance | Student t-test |
| Two groups with unequal variance | Welch t-test |
| Two non-normal groups | Mann-Whitney U |
| Multiple groups | one-way ANOVA |
| Multiple genes and multiple tests | Benjamini-Hochberg FDR |

Run statistical tests on Delta Ct values rather than directly on fold change values.

## Significance Labels

| Label | Threshold |
|---|---:|
| ns | p or FDR >= 0.05 |
| * | < 0.05 |
| ** | < 0.01 |
| *** | < 0.001 |

## Suggested Outputs

- `summary.csv`: fold change, log2FC, p value, and FDR for each gene
- `detail.csv`: Ct, Delta Ct, and Delta Delta Ct for each sample
- `figures/`: bar plots, scatter plots, or volcano plots
- `provenance.json`: parameters, input file hashes, software version, and run time
