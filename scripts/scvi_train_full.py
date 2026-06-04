# Save as ~/Desktop/sea_ad_analysis/scripts/scvi_train_full.py
import os
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import scvi
import scanpy as sc
import scipy.io
import pandas as pd
import numpy as np

os.makedirs(os.path.expanduser(
    "~/sea_ad_analysis/results"), exist_ok=True)

print("Loading FULL dataset...")

# Load full counts matrix
counts = scipy.io.mmread(
    os.path.expanduser(
        "~/Downloads/sea_ad_counts.mtx")).T.tocsr()

# Load full gene list
genes = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_genes_full.csv"))

# Load metadata
metadata = pd.read_csv(
    os.path.expanduser(
        "~/Downloads/sea_ad_metadata.csv"), index_col=0)

print(f"Counts shape: {counts.shape}")
print(f"Genes: {len(genes)}")
print(f"Metadata: {metadata.shape}")

# Align dimensions
n_cells = min(counts.shape[0], len(metadata))
counts = counts[:n_cells]
metadata = metadata.iloc[:n_cells]

# Create AnnData
adata = sc.AnnData(X=counts)
adata.var_names = genes['gene'].values
adata.obs = metadata.copy()
adata.obs_names = metadata.index

print(f"\nAnnData: {adata.n_obs} cells x {adata.n_vars} genes")

# ── QC filtering ──────────────────────────────────────────────────────────────
print("\nQC filtering...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
print(f"After filtering: {adata.n_obs} cells x {adata.n_vars} genes")

# Save filtered metadata
adata.obs.to_csv(os.path.expanduser(
    "~/Downloads/sea_ad_metadata_full_filtered.csv"))
print("Filtered metadata saved!")

# ── Highly variable genes ─────────────────────────────────────────────────────
# Store raw counts first
adata.layers["counts"] = adata.X.copy()

# Normalize for HVG selection
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Select highly variable genes — use more for full dataset
print("\nSelecting highly variable genes...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,       # more genes than before
    subset=True,            # keep only HVGs
    layer="counts",
    flavor="seurat_v3"
)
print(f"After HVG selection: {adata.n_vars} genes")

# ── Setup scVI ────────────────────────────────────────────────────────────────
print("\nSetting up scVI...")
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["donor_id"],
    continuous_covariate_keys=["Fraction.mitochrondrial.UMIs"]
)

# ── Train scVI ────────────────────────────────────────────────────────────────
print("Training scVI model (this will take longer)...")
model = scvi.model.SCVI(
    adata,
    n_layers=2,
    n_latent=30,        # more latent dims for larger dataset
    gene_likelihood="nb"
)

model.train(
    max_epochs=150,     # more epochs for larger dataset
    early_stopping=True,
    batch_size=512      # larger batch for efficiency
)
print("Training complete!")

# ── Save everything ───────────────────────────────────────────────────────────
print("\nSaving results...")

# Save latent
latent = model.get_latent_representation()
np.save(os.path.expanduser(
    "~/Downloads/scvi_latent_full.npy"), latent)
print(f"Latent saved: {latent.shape}")

# Save model
model.save(os.path.expanduser(
    "~/sea_ad_analysis/results/scvi_model_full"),
    overwrite=True)

# Save gene names used
pd.DataFrame({"gene": adata.var_names}).to_csv(
    os.path.expanduser(
        "~/Downloads/sea_ad_genes_hvg_full.csv"),
    index=False)

# Save updated metadata
adata.obs.to_csv(os.path.expanduser(
    "~/Downloads/sea_ad_metadata_full_filtered.csv"))

print("\n=== Done! ===")
print(f"Cells: {adata.n_obs}")
print(f"HVGs: {adata.n_vars}")
print(f"Latent dims: {latent.shape[1]}")
print("\nNext: run scvi_umap_full.py")
