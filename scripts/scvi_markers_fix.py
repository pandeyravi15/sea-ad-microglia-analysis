# Save as ~/Desktop/sea_ad_analysis/scripts/scvi_markers_fix.py
import os
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io

print("Loading data...")

latent = np.load(os.path.expanduser("~/Downloads/scvi_latent.npy"))
embedding = np.load(os.path.expanduser("~/Downloads/scvi_umap.npy"))
metadata = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata_clustered.csv"),
    index_col=0)

raw_counts = scipy.io.mmread(
    os.path.expanduser("~/Downloads/sea_ad_counts_subset.mtx")).T.tocsr()
genes = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_genes_subset.csv"))

n_cells = len(metadata)
raw_counts = raw_counts[:n_cells]

adata = anndata.AnnData(X=raw_counts)
adata.var_names = genes['gene'].values
adata.obs = metadata.copy()
adata.obsm["X_scVI"] = latent
adata.obsm["X_umap"] = embedding

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Make sure clusters are correct type
adata.obs["scVI_clusters"] = adata.obs["scVI_clusters"].astype(str).astype("category")

print(f"AnnData: {adata.n_obs} cells, {adata.n_vars} genes")
print(f"Clusters: {adata.obs['scVI_clusters'].value_counts().sort_index()}")

# ── Debug step — check what rank_genes returns ────────────────────────────────
print("\nRunning rank_genes_groups...")
sc.tl.rank_genes_groups(
    adata,
    groupby="scVI_clusters",
    method="wilcoxon",
    key_added="markers"
)

# Check raw structure
print("\nKeys in markers result:")
print(list(adata.uns["markers"].keys()))

# Extract markers manually — bypass the DataFrame issue
print("\n=== Top markers per cluster ===")
os.makedirs(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/markers"), exist_ok=True)

result = adata.uns["markers"]
groups = result["names"].dtype.names
print(f"Groups found: {groups}")

all_markers = {}
for cluster in groups:
    names = result["names"][cluster]
    scores = result["scores"][cluster]
    pvals = result["pvals_adj"][cluster]
    logfc = result["logfoldchanges"][cluster]

    df = pd.DataFrame({
        "gene": names,
        "logfoldchange": logfc,
        "pval_adj": pvals,
        "score": scores
    })

    # Sort by logfoldchange
    df = df.sort_values("logfoldchange", ascending=False)
    all_markers[cluster] = df

    print(f"\nCluster {cluster} top 10:")
    print(df.head(10)[["gene","logfoldchange","pval_adj"]].to_string())

    df.to_csv(os.path.expanduser(
        f"~/Desktop/sea_ad_analysis/results/markers/cluster_{cluster}_markers.csv"),
        index=False)

print("\nMarkers saved!")

# ── AD gene expression ────────────────────────────────────────────────────────
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
print(f"\nAD genes present: {list(present.keys())}")

# ── Mean expression per cluster ───────────────────────────────────────────────
print("\n=== Mean expression of AD genes per cluster ===")
expr_rows = []
for cluster in sorted(adata.obs["scVI_clusters"].unique(),
                       key=lambda x: int(x)):
    mask = (adata.obs["scVI_clusters"] == cluster).values
    row = {"cluster": cluster}
    for name, gene_id in present.items():
        if gene_id in adata.var_names:
            idx = list(adata.var_names).index(gene_id)
            vals = np.asarray(adata.X[mask, idx].todense()).flatten()
            row[name] = round(float(vals.mean()), 3)
    expr_rows.append(row)

expr_df = pd.DataFrame(expr_rows).set_index("cluster")
print(expr_df)
expr_df.to_csv(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/AD_gene_mean_expression.csv"))

# ── Heatmap of AD gene expression ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(expr_df.values, aspect="auto", cmap="RdYlBu_r")
ax.set_xticks(range(len(expr_df.columns)))
ax.set_xticklabels(expr_df.columns, rotation=45, ha="right", fontsize=11)
ax.set_yticks(range(len(expr_df.index)))
ax.set_yticklabels([f"Cluster {i}" for i in expr_df.index], fontsize=11)
plt.colorbar(im, ax=ax, label="Mean expression (log1p)")

# Annotate values
for i in range(len(expr_df.index)):
    for j in range(len(expr_df.columns)):
        ax.text(j, i, f"{expr_df.values[i,j]:.2f}",
                ha="center", va="center", fontsize=9,
                color="black")

ax.set_title("AD Gene Expression Across scVI Clusters\n(SEA-AD Human Microglia DLPFC)",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/scVI_AD_heatmap.png"),
    dpi=300, bbox_inches="tight")
plt.close()
print("Heatmap saved!")

# ── Violin plots ──────────────────────────────────────────────────────────────
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
        axes[i].set_xticklabels(
            [f"C{c}" for c in sorted(
                adata.obs["scVI_clusters"].unique(),
                key=lambda x: int(x))])
    plt.suptitle("DAM Gene Expression Across scVI Clusters",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.expanduser(
        "~/Desktop/sea_ad_analysis/results/scVI_DAM_violin.png"),
        dpi=300, bbox_inches="tight")
    plt.close()
    print("Violin plot saved!")

print("\n=== Complete! ===")
print("Files saved:")
print("  ~/Desktop/sea_ad_analysis/results/markers/cluster_*_markers.csv")
print("  ~/Desktop/sea_ad_analysis/results/AD_gene_mean_expression.csv")
print("  ~/Desktop/sea_ad_analysis/results/scVI_AD_heatmap.png")
print("  ~/Desktop/sea_ad_analysis/results/scVI_DAM_violin.png")
