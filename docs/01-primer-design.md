# 01. Primer Design

## Goal

Design specific, stable, and efficient qPCR primers to reduce genomic DNA amplification, primer dimers, nonspecific products, and Ct-value bias.

## Recommended Workflow

```mermaid
flowchart LR
    A["Define target gene"] --> B["Check transcripts and exon structure"]
    B --> C["Initial design with PrimerBank or Primer-BLAST"]
    C --> D["Check Tm, GC content, primer length, and amplicon length"]
    D --> E["Validate specificity with BLAST"]
    E --> F["Order primers"]
    F --> G["Validate amplification efficiency with a standard curve"]
```

## Design Parameters

| Parameter | Recommended range | Notes |
|---|---:|---|
| Primer length | 18-24 bp | Around 20 bp is commonly used. |
| Tm | 58-62 degC | Forward and reverse primers should differ by less than 2 degC. |
| GC content | 40%-60% | Avoid very high GC content, which can impair annealing. |
| Amplicon length | 80-200 bp | 100-150 bp is usually robust for qPCR. |
| 3' end | G/C is acceptable, but avoid long G/C runs | Reduces nonspecific priming risk. |
| Self-complementarity | < 4 bp | Helps avoid hairpin structures. |
| Primer dimer risk | Minimize 3' complementarity | Dimers affect SYBR Green signal. |

## Genomic DNA Control

Preferred strategies:

1. Design primers across an exon-exon junction.
2. Place forward and reverse primers in different exons.
3. Span a long intron so the genomic DNA amplicon is too large to amplify efficiently.
4. Include an NRT control to check genomic DNA carryover.

## Tools

| Tool | Purpose |
|---|---|
| Primer-BLAST | Primer design and specificity checking |
| PrimerBank | Validated human and mouse primer lookup |
| Primer3 | Custom primer design |
| IDT OligoAnalyzer | Hairpin, dimer, and Tm evaluation |

## Pre-Run Checklist

- [ ] Target gene and species are correct.
- [ ] Transcript accession is recorded.
- [ ] Primers span exons or can distinguish cDNA from genomic DNA.
- [ ] Primer-BLAST shows no major off-target products.
- [ ] Amplicon length is suitable for qPCR.
- [ ] Standard-curve efficiency is within 90%-110%.
