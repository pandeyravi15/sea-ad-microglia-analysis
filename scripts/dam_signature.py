# Save as ~/Desktop/sea_ad_analysis/scripts/dam_signature.py
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

adata.obs["scVI_clusters"] = adata.obs["scVI_clusters"].astype(str).astype("category")

# ── Comprehensive AD gene signatures ─────────────────────────────────────────
# DAM signature genes (Keren-Shaul et al. 2017 + human validation)
dam_genes = {
    # Homeostatic markers — should DECREASE in DAM
    "Homeostatic": {
        "P2RY12":  "ENSG00000169313",
        "CX3CR1":  "ENSG00000168329",
        "TMEM119": "ENSG00000183160",
        "CSF1R":   "ENSG00000182578",
    },
    # DAM markers — should INCREASE in disease
    "DAM_activated": {
        "TREM2":  "ENSG00000095970",
        "TYROBP": "ENSG00000011600",
        "SPP1":   "ENSG00000118785",
        "CD44":   "ENSG00000026508",
        "CR1":    "ENSG00000203710",
        "LPL":    "ENSG00000175445",
        "CST7":   "ENSG00000077984",
        "APOE":   "ENSG00000130203",
    }
}

# Filter to present genes
print("\nGenes present in dataset:")
for group, genes_dict in dam_genes.items():
    present = {k: v for k, v in genes_dict.items()
               if v in adata.var_names}
    dam_genes[group] = present
    print(f"  {group}: {list(present.keys())}")

# ── Compute DAM score per cell ────────────────────────────────────────────────
print("\nComputing DAM scores...")

# DAM activation score
dam_gene_ids = list(dam_genes["DAM_activated"].values())
homeostatic_gene_ids = list(dam_genes["Homeostatic"].values())

# Score = mean expression of DAM genes - mean expression of homeostatic genes
dam_expr = np.zeros(adata.n_obs)
for gene_id in dam_gene_ids:
    if gene_id in adata.var_names:
        idx = list(adata.var_names).index(gene_id)
        dam_expr += np.asarray(
            adata.X[:, idx].todense()).flatten()
dam_expr /= max(len(dam_gene_ids), 1)

homeostatic_expr = np.zeros(adata.n_obs)
for gene_id in homeostatic_gene_ids:
    if gene_id in adata.var_names:
        idx = list(adata.var_names).index(gene_id)
        homeostatic_expr += np.asarray(
            adata.X[:, idx].todense()).flatten()
homeostatic_expr /= max(len(homeostatic_gene_ids), 1)

# DAM score = activation - homeostatic
adata.obs["DAM_score"] = dam_expr - homeostatic_expr
adata.obs["DAM_activation"] = dam_expr
adata.obs["Homeostatic_score"] = homeostatic_expr

print(f"\nDAM score per cluster:")
score_summary = adata.obs.groupby("scVI_clusters")[
    ["DAM_score","DAM_activation","Homeostatic_score"]].mean()
print(score_summary.round(3))
score_summary.to_csv(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/DAM_scores_per_cluster.csv"))

# ── Plot 1 — DAM score on UMAP ───────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# DAM score
sc_plot = axes[0].scatter(
    embedding[:, 0], embedding[:, 1],
    c=adata.obs["DAM_score"].values,
    cmap="RdYlBu_r", s=1, alpha=0.5)
plt.colorbar(sc_plot, ax=axes[0], label="DAM score")
axes[0].set_title("DAM Score\n(Activation - Homeostatic)",
                   fontsize=12, fontweight="bold")
axes[0].set_xlabel("UMAP 1")
axes[0].set_ylabel("UMAP 2")

# DAM activation
sc_plot2 = axes[1].scatter(
    embedding[:, 0], embedding[:, 1],
    c=adata.obs["DAM_activation"].values,
    cmap="Reds", s=1, alpha=0.5)
plt.colorbar(sc_plot2, ax=axes[1], label="Expression")
axes[1].set_title("DAM Activation Score\n(TREM2/TYROBP/SPP1/CD44)",
                   fontsize=12, fontweight="bold")
axes[1].set_xlabel("UMAP 1")

# Homeostatic
sc_plot3 = axes[2].scatter(
    embedding[:, 0], embedding[:, 1],
    c=adata.obs["Homeostatic_score"].values,
    cmap="Blues", s=1, alpha=0.5)
plt.colorbar(sc_plot3, ax=axes[2], label="Expression")
axes[2].set_title("Homeostatic Score\n(P2RY12/CX3CR1/TMEM119/CSF1R)",
                   fontsize=12, fontweight="bold")
axes[2].set_xlabel("UMAP 1")

plt.suptitle("DAM vs Homeostatic Microglia Signatures\nSEA-AD Human DLPFC",
             fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/DAM_score_UMAP.png"),
    dpi=300, bbox_inches="tight")
plt.close()
print("DAM score UMAP saved!")

# ── Plot 2 — DAM score by disease group ──────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

disease_order = ["Reference", "Control", "AD_Low_Int", "AD_High"]
colors = ["grey", "steelblue", "orange", "red"]

# Violin — DAM score by disease
data_by_disease = [
    adata.obs[adata.obs["disease_group"] == g]["DAM_score"].values
    for g in disease_order
    if g in adata.obs["disease_group"].values
]
labels = [g for g in disease_order
          if g in adata.obs["disease_group"].values]
cols = [c for g, c in zip(disease_order, colors)
        if g in adata.obs["disease_group"].values]

vp = axes[0].violinplot(data_by_disease, positions=range(len(labels)),
                         showmedians=True)
for patch, color in zip(vp['bodies'], cols):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
axes[0].set_xticks(range(len(labels)))
axes[0].set_xticklabels(labels, rotation=45, ha="right")
axes[0].set_title("DAM Score by Disease Group", fontsize=12,
                   fontweight="bold")
axes[0].set_ylabel("DAM Score")
axes[0].axhline(y=0, color="black", linestyle="--", alpha=0.5)

# Violin — DAM score by cluster
cluster_order = sorted(adata.obs["scVI_clusters"].unique(),
                        key=lambda x: int(x))
data_by_cluster = [
    adata.obs[adata.obs["scVI_clusters"] == c]["DAM_score"].values
    for c in cluster_order
]
vp2 = axes[1].violinplot(data_by_cluster,
                          positions=range(len(cluster_order)),
                          showmedians=True)
axes[1].set_xticks(range(len(cluster_order)))
axes[1].set_xticklabels([f"Cluster {c}" for c in cluster_order],
                         rotation=45, ha="right")
axes[1].set_title("DAM Score by scVI Cluster", fontsize=12,
                   fontweight="bold")
axes[1].set_ylabel("DAM Score")
axes[1].axhline(y=0, color="black", linestyle="--", alpha=0.5)

plt.suptitle("DAM Activation Across Disease Groups and Clusters",
             fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/DAM_score_violin.png"),
    dpi=300, bbox_inches="tight")
plt.close()
print("Violin plot saved!")

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n=== DAM Score by Disease Group ===")
print(adata.obs.groupby("disease_group")["DAM_score"].mean().round(3))

print("\n=== DAM Score by Cluster ===")
print(adata.obs.groupby("scVI_clusters")["DAM_score"].mean().round(3))

print("\nDone! Files saved to ~/Desktop/sea_ad_analysis/results/")
