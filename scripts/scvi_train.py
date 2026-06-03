import scvi
import scanpy as sc
import scipy.io
import pandas as pd
import numpy as np
import os

os.makedirs(os.path.expanduser("~/Desktop/sea_ad_analysis/results"), exist_ok=True)

print("Loading data...")
counts = scipy.io.mmread(
    os.path.expanduser("~/Downloads/sea_ad_counts_subset.mtx")).T.tocsr()
genes = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_genes_subset.csv"))
metadata = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata.csv"), index_col=0)

adata = sc.AnnData(X=counts)
adata.var_names = genes['gene'].values
adata.obs = metadata
adata.obs_names = metadata.index

# Filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")

# Save filtered metadata immediately
adata.obs.to_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata_filtered.csv"))
print("Filtered metadata saved!")

# Store raw counts
adata.layers["counts"] = adata.X.copy()

# Normalize for later
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Setup and train scVI
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["donor_id"],
    continuous_covariate_keys=["Fraction.mitochrondrial.UMIs"]
)

model = scvi.model.SCVI(adata, n_layers=2, n_latent=20, gene_likelihood="nb")
model.train(max_epochs=100, early_stopping=True, batch_size=256)
print("Training complete!")

# Save latent representation
latent = model.get_latent_representation()
np.save(os.path.expanduser("~/Downloads/scvi_latent.npy"), latent)
print(f"Latent saved! Shape: {latent.shape}")

# Save model
model.save(
    os.path.expanduser("~/Desktop/sea_ad_analysis/results/scvi_model"),
    overwrite=True)
print("Model saved!")
print("\nNow run: python3 scvi_umap.py")
