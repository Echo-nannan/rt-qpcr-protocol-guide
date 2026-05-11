# 03. Reverse Transcription

## Goal

Reverse-transcribe equal amounts of RNA into cDNA so samples remain comparable during qPCR.

## Calculation Rule

Use the same RNA input amount across one experimental batch, such as 500 ng or 1000 ng per sample.

```text
RNA volume (uL) = target RNA input (ng) / RNA concentration (ng/uL)
```

If the calculated RNA volume is too large, concentrate RNA or lower the unified RNA input for all samples.

## Example 20 uL RT Reaction

| Component | Volume | Notes |
|---|---:|---|
| RNA template | X uL | Equal total RNA input |
| 4x gDNA wiper Mix | 4 uL | Genomic DNA removal |
| RNase-free water | To 16 uL | First-stage reaction volume |
| 5x qRT SuperMix | 4 uL | Reverse-transcription mix |
| Total volume | 20 uL |  |

## Example Program

| Step | Temperature | Time |
|---|---:|---:|
| gDNA removal | 42 degC | 2 min |
| Primer annealing | 25 degC | 5 min |
| Reverse transcription | 50 degC | 15 min |
| Enzyme inactivation | 85 degC | 5 sec |
| Hold | 4 degC | hold |

## cDNA Handling

- Store cDNA at 4 degC for short-term use and at -20 degC for long-term storage.
- Use a consistent dilution, commonly 5x or 10x, before qPCR.
- Apply the same dilution factor to all samples in one experiment.
- Avoid repeated freeze-thaw cycles; aliquot if needed.
