# Save as ~/Desktop/sea_ad_analysis/scripts/scvi_cluster.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import leidenalg

print("Loading scVI results...")

# Load saved results
latent = np.load(os.path.expanduser("~/Downloads/scvi_latent.npy"))
embedding = np.load(os.path.expanduser("~/Downloads/scvi_umap.npy"))
metadata = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata_filtered.csv"),
    index_col=0)

print(f"Cells: {len(metadata)}")

# Clustering using leiden on scVI latent space
print("\nClustering on scVI latent space...")
import scanpy as sc
import anndata

# Create minimal AnnData for clustering only
adata = anndata.AnnData(X=latent)
adata.obsm["X_scVI"] = latent
adata.obsm["X_umap"] = embedding
adata.obs = metadata.copy()

# Compute neighbors on scVI latent
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)

# Leiden clustering
sc.tl.leiden(adata, resolution=0.3, key_added="scVI_clusters")

print(f"\nClusters found:\n{adata.obs['scVI_clusters'].value_counts().sort_index()}")

# Save updated metadata with clusters
adata.obs.to_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata_clustered.csv"))
print("Clustered metadata saved!")

# Plot clusters on UMAP
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# Get unique clusters and colors
clusters = adata.obs["scVI_clusters"].unique()
n_clusters = len(clusters)
colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))
color_map = dict(zip(sorted(clusters, key=lambda x: int(x)), colors))

# Plot 1 — scVI clusters
for cluster, color in color_map.items():
    mask = (adata.obs["scVI_clusters"] == cluster).values
    axes[0].scatter(embedding[mask, 0], embedding[mask, 1],
                   c=[color], label=cluster, s=1, alpha=0.5)
axes[0].legend(markerscale=5, fontsize=8,
               title="Cluster", bbox_to_anchor=(1.05, 1))
axes[0].set_title("scVI Clusters", fontsize=13)
axes[0].set_xlabel("UMAP 1")
axes[0].set_ylabel("UMAP 2")

# Plot 2 — Disease group overlay
disease_colors = {
    "Reference":  "grey",
    "Control":    "steelblue",
    "AD_Low_Int": "orange",
    "AD_High":    "red"
}
for group, color in disease_colors.items():
    mask = (metadata["disease_group"] == group).values
    axes[1].scatter(embedding[mask, 0], embedding[mask, 1],
                   c=color, label=group, s=1, alpha=0.5)
axes[1].legend(markerscale=5, fontsize=9)
axes[1].set_title("Disease Group", fontsize=13)
axes[1].set_xlabel("UMAP 1")

# Plot 3 — Cluster x Disease heatmap
cluster_disease = pd.crosstab(
    adata.obs["scVI_clusters"],
    adata.obs["disease_group"],
    normalize="index"
) * 100

im = axes[2].imshow(cluster_disease.values, aspect="auto", cmap="RdYlBu_r")
axes[2].set_xticks(range(len(cluster_disease.columns)))
axes[2].set_xticklabels(cluster_disease.columns, rotation=45, ha="right")
axes[2].set_yticks(range(len(cluster_disease.index)))
axes[2].set_yticklabels([f"Cluster {i}" for i in cluster_disease.index])
axes[2].set_title("% Disease per Cluster", fontsize=13)
plt.colorbar(im, ax=axes[2], label="% cells")

# Annotate heatmap
for i in range(len(cluster_disease.index)):
    for j in range(len(cluster_disease.columns)):
        axes[2].text(j, i, f"{cluster_disease.values[i,j]:.0f}",
                    ha="center", va="center", fontsize=8)

plt.suptitle("SEA-AD Microglia — scVI Clustering\nRavi Shanker Pandey, Ph.D.",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(
    os.path.expanduser("~/Desktop/sea_ad_analysis/results/scVI_clusters.png"),
    dpi=300, bbox_inches="tight")
plt.close()
print("Cluster plot saved!")

# Identify DAM clusters — high AD_High proportion
print("\n=== DAM Cluster Identification ===")
print("Clusters enriched for AD_High (>50%):")
if "AD_High" in cluster_disease.columns:
    dam_clusters = cluster_disease[
        cluster_disease["AD_High"] > 50].index.tolist()
    print(f"DAM clusters: {dam_clusters}")
    
print("\nFull cluster-disease table:")
print(cluster_disease.round(1))

# Save cluster summary
cluster_disease.to_csv(
    os.path.expanduser(
        "~/Desktop/sea_ad_analysis/results/scVI_cluster_disease_pct.csv"))
print("\nDone!")
