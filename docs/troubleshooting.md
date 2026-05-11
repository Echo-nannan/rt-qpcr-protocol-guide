# Troubleshooting

| Problem | Common causes | Suggested actions |
|---|---|---|
| No amplification | Missing template, wrong primers, failed mix | Check cDNA, primer sequence, reaction setup, and cycling program. |
| Ct > 35 | Low template input, low target expression, poor efficiency | Increase template input, reduce dilution, or optimize primers. |
| Ct is too low | Too much template, contamination, or incorrect threshold | Dilute template and check NTC and threshold settings. |
| Technical replicate SD > 0.5 | Pipetting error, bubbles, poor mixing, sealing issue | Spin plate again and review well positions and pipetting records. |
| NTC amplification | Reagent contamination or primer dimers | Replace water and mix, remake primers, and inspect melt curves. |
| NRT amplification | Genomic DNA carryover | Strengthen DNase or gDNA-wiper treatment and redesign exon-spanning primers. |
| Multiple melt-curve peaks | Nonspecific amplification or primer dimers | Raise annealing temperature, lower primer concentration, or redesign primers. |
| Efficiency < 90% | Poor primer efficiency, inhibitors, or low template quality | Run a standard curve, dilute template, or re-extract RNA. |
| Efficiency > 110% | Primer dimers or standard-curve error | Inspect melt curves and repeat serial dilution. |
| Unstable reference gene | Reference gene affected by treatment or uneven sample quality | Replace or add reference genes and use multi-reference normalization. |

## Troubleshooting Order

1. Check NTC and NRT first to rule out contamination and genomic DNA.
2. Review melt curves to rule out nonspecific amplification.
3. Review technical replicate SD to catch pipetting or well-position issues.
4. Check reference gene stability to evaluate sample quality and normalization reliability.
5. Review statistical settings and Delta Delta Ct parameters last.
