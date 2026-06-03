# Save as ~/sea_ad_analysis/scripts/scvi_markers.py

import os
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io

print("Loading data for marker analysis...")

# Load all saved results
latent = np.load(os.path.expanduser("~/Downloads/scvi_latent.npy"))
embedding = np.load(os.path.expanduser("~/Downloads/scvi_umap.npy"))
metadata = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata_clustered.csv"),
    index_col=0)

# Load normalized counts
raw_counts = scipy.io.mmread(
    os.path.expanduser("~/Downloads/sea_ad_counts_subset.mtx")).T.tocsr()
genes = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_genes_subset.csv"))

# Align cells
n_cells = len(metadata)
raw_counts = raw_counts[:n_cells]

# Create AnnData
adata = anndata.AnnData(X=raw_counts)
adata.var_names = genes['gene'].values
adata.obs = metadata.copy()
adata.obsm["X_scVI"] = latent
adata.obsm["X_umap"] = embedding

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Set clusters
adata.obs["scVI_clusters"] = adata.obs["scVI_clusters"].astype("category")

print(f"AnnData: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Clusters: {adata.obs['scVI_clusters'].value_counts().sort_index()}")


# ── Fix 1: Marker genes without pval filter first ────────────────────────────
print("\nFinding markers...")
sc.tl.rank_genes_groups(adata,
                         groupby="scVI_clusters",
                         method="wilcoxon",
                         key_added="markers",
                         pts=True)

# Check raw results first
print("\n=== Top markers per cluster (no filter) ===")
os.makedirs(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/markers"), exist_ok=True)

for cluster in sorted(adata.obs["scVI_clusters"].unique(), 
                       key=lambda x: int(x)):
    df = sc.get.rank_genes_groups_df(
        adata, group=cluster, key="markers")
    print(f"\nCluster {cluster} top 5:")
    print(df[["names","logfoldchanges","pvals_adj"]].head(5).to_string())
    df.to_csv(os.path.expanduser(
        f"~/Desktop/sea_ad_analysis/results/markers/cluster_{cluster}_markers.csv"),
        index=False)


# ── AD gene expression across clusters ───────────────────────────────────────
ad_genes = {
    "TREM2":   "ENSG00000095970",
    "TYROBP":  "ENSG00000011600",
    "SPP1":    "ENSG00000118785",
    "P2RY12":  "ENSG00000169313",
    "CX3CR1":  "ENSG00000168329",
    "CSF1R":   "ENSG00000182578",
    "TMEM119": "ENSG00000183160",
    "CD44":    "ENSG00000026508",
    "CR1":     "ENSG00000203710"
}

# ── Fix 2: AD genes dotplot ──────────────────────────────────────────────────
ad_genes = {
    "TREM2":   "ENSG00000095970",
    "TYROBP":  "ENSG00000011600",
    "SPP1":    "ENSG00000118785",
    "P2RY12":  "ENSG00000169313",
    "CX3CR1":  "ENSG00000168329",
    "CSF1R":   "ENSG00000182578",
    "TMEM119": "ENSG00000183160",
    "CD44":    "ENSG00000026508",
    "CR1":     "ENSG00000203710"
}

present = {k: v for k, v in ad_genes.items() if v in adata.var_names}
print(f"\nAD genes found: {list(present.keys())}")

# Dotplot — save directly without axis manipulation
dp = sc.pl.dotplot(adata,
                    var_names=list(present.values()),
                    groupby="scVI_clusters",
                    show=False,
                    return_fig=True,
                    title="AD Gene Expression Across scVI Clusters",
                    var_group_labels=list(present.keys()),
                    var_group_positions=[(i, i) for i in
                                         range(len(present))])
dp.savefig(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/scVI_AD_dotplot.png"),
    dpi=300, bbox_inches="tight")
plt.close()
print("Dotplot saved!")

# ── Violin plots ─────────────────────────────────────────────────────────────
plot_genes = {k: v for k, v in present.items()
              if k in ["TYROBP", "SPP1", "P2RY12", "CX3CR1"]}

if plot_genes:
    fig, axes = plt.subplots(1, len(plot_genes), figsize=(16, 5))
    if len(plot_genes) == 1:
        axes = [axes]
    for i, (name, gene_id) in enumerate(plot_genes.items()):
        sc.pl.violin(adata,
                     keys=gene_id,
                     groupby="scVI_clusters",
                     ax=axes[i],
                     show=False,
                     rotation=0)
        axes[i].set_title(name, fontsize=13, fontweight="bold")
        axes[i].set_xlabel("Cluster")
    plt.suptitle("DAM Gene Expression — scVI Clusters",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.expanduser(
        "~/Desktop/sea_ad_analysis/results/scVI_DAM_violin.png"),
        dpi=300, bbox_inches="tight")
    plt.close()
    print("Violin plot saved!")
    
    
print("\n=== Analysis Complete ===")

# ── Manual mean expression per cluster ──────────────────────────────────────
print("\n=== Mean expression of AD genes per cluster ===")
expr_df = pd.DataFrame(index=adata.obs_names)
for name, gene_id in present.items():
    if gene_id in adata.var_names:
        idx = list(adata.var_names).index(gene_id)
        expr_df[name] = np.asarray(adata.X[:, idx].todense()).flatten()

expr_df["cluster"] = adata.obs["scVI_clusters"].values
mean_expr = expr_df.groupby("cluster").mean()
print(mean_expr.round(3))
mean_expr.to_csv(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/AD_gene_mean_expression.csv"))

print("\nDone!")
