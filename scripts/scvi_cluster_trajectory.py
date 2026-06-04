# Save as ~/Desktop/sea_ad_analysis/scripts/scvi_cluster_full.py
import os
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata

print("Loading full dataset...")

latent = np.load(os.path.expanduser("~/Downloads/scvi_latent_full.npy"))
embedding = np.load(os.path.expanduser("~/Downloads/scvi_umap_full.npy"))
metadata = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata_full_clustered.csv"),
    index_col=0)

print(f"Cells: {len(metadata)}")

# Create AnnData
adata = anndata.AnnData(X=latent)
adata.obsm["X_scVI"] = latent
adata.obsm["X_umap"] = embedding
adata.obs = metadata.copy()

adata.obs["scVI_clusters_full"] = (adata.obs["scVI_clusters_full"]
                                    .astype(str)
                                    .astype("category"))

# ── Plot clusters + disease + heatmap ────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(22, 7))

# Cluster colors
cluster_colors = {
    "0": "red",       "1": "orange",
    "2": "steelblue", "3": "green",
    "4": "purple",    "5": "brown"
}

# Plot 1 — clusters
for cluster, color in cluster_colors.items():
    mask = (adata.obs["scVI_clusters_full"] == cluster).values
    n = mask.sum()
    axes[0].scatter(embedding[mask,0], embedding[mask,1],
                   c=color, label=f"C{cluster} (n={n})",
                   s=0.5, alpha=0.5)
axes[0].legend(markerscale=6, fontsize=9)
axes[0].set_title("scVI Clusters — Full Dataset",
                   fontsize=13, fontweight="bold")
axes[0].set_xlabel("UMAP 1")
axes[0].set_ylabel("UMAP 2")

# Plot 2 — disease
disease_colors = {
    "Reference":  "grey",   "Control":    "steelblue",
    "AD_Low_Int": "orange", "AD_High":    "red"
}
for group, color in disease_colors.items():
    mask = (metadata["disease_group"] == group).values
    axes[1].scatter(embedding[mask,0], embedding[mask,1],
                   c=color, label=group, s=0.5, alpha=0.4)
axes[1].legend(markerscale=6, fontsize=9)
axes[1].set_title("Disease Group", fontsize=13, fontweight="bold")
axes[1].set_xlabel("UMAP 1")

# Plot 3 — heatmap
cluster_disease = pd.crosstab(
    adata.obs["scVI_clusters_full"],
    adata.obs["disease_group"],
    normalize="index") * 100

im = axes[2].imshow(cluster_disease.values,
                     aspect="auto", cmap="RdYlBu_r")
axes[2].set_xticks(range(len(cluster_disease.columns)))
axes[2].set_xticklabels(cluster_disease.columns,
                         rotation=45, ha="right")
axes[2].set_yticks(range(len(cluster_disease.index)))
axes[2].set_yticklabels([f"Cluster {i}"
                          for i in cluster_disease.index])
axes[2].set_title("% Disease per Cluster",
                   fontsize=13, fontweight="bold")
plt.colorbar(im, ax=axes[2], label="% cells")

for i in range(len(cluster_disease.index)):
    for j in range(len(cluster_disease.columns)):
        axes[2].text(j, i,
                     f"{cluster_disease.values[i,j]:.0f}",
                     ha="center", va="center", fontsize=9)

plt.suptitle("SEA-AD Microglia — Full Dataset scVI Clustering\n"
             "Ravi Shanker Pandey, Ph.D.",
             fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/scVI_clusters_full.png"),
    dpi=300, bbox_inches="tight")
plt.close()
print("Cluster plot saved!")

# ── DAM signature on full dataset ────────────────────────────────────────────
import scipy.io

print("\nLoading counts for DAM scoring...")
raw_counts = scipy.io.mmread(
    os.path.expanduser(
        "~/Downloads/sea_ad_counts_subset.mtx")).T.tocsr()
genes = pd.read_csv(
    os.path.expanduser(
        "~/Downloads/sea_ad_genes_subset.csv"))

# Use subset counts aligned to full metadata
n_cells = min(raw_counts.shape[0], len(metadata))
raw_counts = raw_counts[:n_cells]
metadata_sub = metadata.iloc[:n_cells]

adata_expr = anndata.AnnData(X=raw_counts)
adata_expr.var_names = genes['gene'].values
adata_expr.obs = metadata_sub.copy()

sc.pp.normalize_total(adata_expr, target_sum=1e4)
sc.pp.log1p(adata_expr)

# DAM genes
dam_genes = {
    "TREM2":   "ENSG00000095970",
    "TYROBP":  "ENSG00000011600",
    "SPP1":    "ENSG00000118785",
    "P2RY12":  "ENSG00000169313",
    "CX3CR1":  "ENSG00000168329",
    "CSF1R":   "ENSG00000182578"
}
present = {k: v for k, v in dam_genes.items()
           if v in adata_expr.var_names}
print(f"DAM genes found: {list(present.keys())}")

# Score
dam_expr = np.zeros(adata_expr.n_obs)
homeo_expr = np.zeros(adata_expr.n_obs)
dam_ids = [v for k, v in present.items()
           if k in ["TREM2","TYROBP","SPP1"]]
homeo_ids = [v for k, v in present.items()
             if k in ["P2RY12","CX3CR1","CSF1R"]]

for g in dam_ids:
    idx = list(adata_expr.var_names).index(g)
    dam_expr += np.asarray(
        adata_expr.X[:,idx].todense()).flatten()
dam_expr /= max(len(dam_ids), 1)

for g in homeo_ids:
    idx = list(adata_expr.var_names).index(g)
    homeo_expr += np.asarray(
        adata_expr.X[:,idx].todense()).flatten()
homeo_expr /= max(len(homeo_ids), 1)

adata_expr.obs["DAM_score"] = dam_expr - homeo_expr

print("\n=== DAM Score by Disease Group ===")
print(adata_expr.obs.groupby(
    "disease_group")["DAM_score"].mean().round(3))

print("\n=== DAM Score by Cluster ===")
print(adata_expr.obs.groupby(
    "scVI_clusters_full")["DAM_score"].mean().round(3))

print("\nDone!")
