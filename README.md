
# SEA-AD Human Microglia Analysis
## Cross-Species Validation of Alzheimer's Disease Microglial Signatures

**Ravi Shanker Pandey, Ph.D.**  
Associate Computational Scientist | The Jackson Laboratory  
[GitHub](https://github.com/pandeyravi15) | [ORCID](https://orcid.org/0000-0001-9567-2851)

---

## Overview

Comprehensive single-cell RNA-seq analysis of 42,486 human microglia from the 
Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD) — validating disease-associated 
microglial (DAM) signatures identified in mouse AD models against human disease.

---

## Key Findings

### 1. scVI Probabilistic Integration
- Trained Variational Autoencoder on 42,486 microglia from **42 donors**
- Explicit donor batch correction via categorical covariates
- 4,000 highly variable genes, 30-dimensional latent space
- Identified **6 distinct microglial clusters** — Cluster 0 most disease-specific (72% AD_High)

### 2. Disease-Associated Microglial States
- **Cluster 0** (n=1,274): 72% AD_High, 0% Reference — primary DAM population
- **Cluster 1** (n=672): 60% AD_High — secondary DAM state
- **Cluster 2** (n=35,544): Large transitional population
- Progressive DAM activation confirmed by SPP1↑, P2RY12↓, CX3CR1↓

### 3. Cross-Species Validation
DAM markers identified in 5XFAD mouse models validated in human SEA-AD data:

| Gene | Mouse (5XFAD) | Human (SEA-AD) |
|------|--------------|----------------|
| TYROBP | ↑ Upregulated | ↑ AD_High enriched |
| SPP1 | ↑ Upregulated | ↑ Progressive increase |
| CSF1R | ↑ Upregulated | ↑ Disease clusters |
| P2RY12 | ↓ Downregulated | ↓ Progressive decrease |
| CX3CR1 | ↓ Downregulated | ↓ Early decrease |

### 4. Slingshot Trajectory Analysis
- Two divergent lineages from common homeostatic origin:
  - **Lineage 1** (4→2→5→0): Primary DAM activation path
  - **Lineage 2** (4→2→3→1): Alternative activation path
- **Novel finding:** Biphasic pseudotime relationship with Braak staging:
  - **Early wave** (Braak 0): High pseudotime — DAM activation precedes tau pathology
  - **Mid-stage** (Braak II-IV): Reduced activation
  - **Late wave** (Braak V-VI): Re-activation with neurodegeneration
- Suggests two therapeutic intervention windows

---

## Repository Structure

```
sea-ad-microglia-analysis/
├── scripts/
│   ├── scvi_train.py           # scVI training — subset (24K cells)
│   ├── scvi_train_full.py      # scVI training — full (42K cells)
│   ├── scvi_umap.py            # UMAP computation
│   ├── scvi_umap_full.py       # UMAP — full dataset
│   ├── scvi_cluster.py         # Leiden clustering
│   ├── scvi_cluster_full.py    # Clustering — full dataset
│   ├── scvi_markers_fix.py     # Marker gene identification
│   ├── dam_signature.py        # DAM signature scoring
│   ├── convert_markers.py      # ENSEMBL to gene symbol conversion
│   └── utils.py                # Shared utility functions
├── plots/
│   ├── SEA_AD_microglia_UMAP.png           # Seurat UMAP
│   ├── SEA_AD_homeostatic_vs_DAM.png       # DAM vs homeostatic genes
│   ├── SEA_AD_disease_trajectory.png       # Gene expression trajectory
│   ├── scVI_UMAP.png                       # scVI UMAP — subset
│   ├── scVI_UMAP_full.png                  # scVI UMAP — full dataset
│   ├── scVI_clusters.png                   # Cluster analysis — subset
│   ├── scVI_clusters_full.png              # Cluster analysis — full
│   ├── DAM_score_UMAP.png                  # DAM scoring on UMAP
│   ├── DAM_score_violin.png                # DAM scores by disease group
│   ├── slingshot_pseudotime.png            # Trajectory — subset
│   ├── slingshot_pseudotime_violin.png     # Pseudotime violins
│   ├── slingshot_full_trajectory.png       # Trajectory — full dataset
│   └── slingshot_full_violin.png           # Full dataset pseudotime
└── results/
├── DAM_scores_per_cluster.csv          # DAM activation scores
├── AD_gene_mean_expression.csv         # AD gene expression summary
├── scVI_cluster_disease_pct.csv        # Disease % per cluster
├── slingshot_pseudotime.csv            # Pseudotime — subset
└── slingshot_full_pseudotime.csv       # Pseudotime — full dataset

```
---

## Methods

### Data Source
- **Dataset:** Microglia-PVM DLPFC — Seattle Alzheimer's Disease Atlas (SEA-AD)
- **Cells:** 42,486 | **Genes:** 35,483 | **Donors:** 42
- **Download:** [CellxGene](https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30)
- **Reference:** Gabitto et al., Nature Neuroscience 2024

### Computational Pipeline

```
Data loading and QC (Seurat/scanpy)
├── Doublet removal
├── min_genes=200, min_cells=10
└── Mitochondrial UMI correction
scVI Integration (scvi-tools v1.4.3)
├── 4,000 highly variable genes
├── Donor batch correction (categorical covariate)
├── Negative binomial gene likelihood
├── n_latent=30, n_layers=2
└── Early stopping, max_epochs=150
Dimensionality reduction
├── UMAP on scVI latent space
└── n_neighbors=15, min_dist=0.3
Clustering (Leiden)
├── KNN graph on scVI latent
├── resolution=0.3
└── igraph implementation
DAM Signature Scoring
├── Keren-Shaul et al. 2017 gene sets
├── DAM: TREM2, TYROBP, SPP1, CD44, CR1
└── Homeostatic: P2RY12, CX3CR1, TMEM119, CSF1R
Trajectory Analysis (Slingshot)
├── Start: Cluster 4 (homeostatic)
├── End: Cluster 0 (primary DAM)
└── Two divergent lineages identified

```
---

### Requirements

```bash
# Python
pip install scvi-tools scanpy umap-learn 
pip install leidenalg igraph anndata scipy

# R
BiocManager::install(c("Seurat", "slingshot",
                        "SingleCellExperiment",
                        "org.Hs.eg.db"))
```

---

## Connection to Mouse Model Work

This analysis directly validates findings from preclinical AD mouse model studies:

- **5XFAD proteomics:** Brain/CSF/plasma signatures validated in human microglia
- **Trem2*R47H splicing:** TREM2 pathway activation confirmed in human DAM clusters  
- **MODEL-AD cross-species framework:** Mouse model signatures recapitulated in SEA-AD

---

## Citation

If you use this analysis, please cite:
- **SEA-AD:** Gabitto et al., Nature Neuroscience 2024
- **scVI:** Lopez et al., Nature Methods 2018
- **Slingshot:** Street et al., BMC Genomics 2018

---

## Related Repositories

| Repository | Description |
|---|---|
| [AD-mouse-splicing-workflow](https://github.com/pandeyravi15/AD-mouse-splicing-workflow) | Cryptic splicing in Trem2*R47H models |
| [Brain-multiomics-pipeline](https://github.com/pandeyravi15/Brain-multiomics-pipeline) | DIA proteomics pipeline |
| [spatial_transcriptomics](https://github.com/pandeyravi15/spatial_transcriptomics) | Human brain Visium analysis |
| [rnaseq-nextflow-aws](https://github.com/pandeyravi15/rnaseq-nextflow-aws) | Cloud-based RNA-seq pipeline |
