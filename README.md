# SEA-AD Human Microglia Analysis

Cross-species validation of Alzheimer's disease microglial signatures 
using the Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD).

## Overview
Analysis of 42,486 microglia from human DLPFC across AD and control 
donors. Validates mouse model AD signatures in human disease.

## Key Findings
- Homeostatic markers (P2RY12, CX3CR1, TMEM119) decrease 
  progressively with AD severity — loss of surveillance microglia
- DAM markers (SPP1, TYROBP, CSF1R) increase progressively 
  with AD neuropathological change
- Cross-species validation confirms 5XFAD mouse model signatures 
  translate to human AD microglia
- Disease-associated microglial clusters (Micro-PVM_3-SEAAD) 
  enriched in high ADNC and late Braak stage donors

## Data Source
Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)
- Dataset: Microglia-PVM DLPFC
- Cells: 42,486 | Genes: 35,483
- Download: https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30
- Reference: Gabitto et al., Nature Neuroscience 2024

## Repository Structure
sea_ad_analysis/
├── scripts/
│   └── analysis.R          # Complete analysis pipeline
├── plots/
│   ├── SEA_AD_microglia_UMAP.png
│   ├── SEA_AD_key_genes_symbols.png
│   ├── SEA_AD_homeostatic_vs_DAM.png
│   └── SEA_AD_disease_trajectory.png
└── results/
├── gene_expression_summary.csv
└── DAM_markers.csv

## Requirements
```r
library(Seurat)
library(zellkonverter)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)
```

## Cross-Species Validation
Key genes validated between 5XFAD mouse model and human SEA-AD:

| Gene | Mouse (5XFAD) | Human (SEA-AD) |
|------|--------------|----------------|
| TYROBP | ↑ Upregulated | ↑ Increases with ADNC |
| SPP1 | ↑ Upregulated | ↑ Progressive increase |
| CSF1R | ↑ Upregulated | ↑ AD_High enriched |
| P2RY12 | ↓ Downregulated | ↓ Progressive decrease |
| CX3CR1 | ↓ Downregulated | ↓ Early decrease |

## Author
Ravi Shanker Pandey, Ph.D.
Associate Computational Scientist, The Jackson Laboratory
[ORCID](https://orcid.org/0000-0001-9567-2851) | 
[GitHub](https://github.com/pandeyravi15)
