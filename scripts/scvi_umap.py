import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
import os

print("Loading saved results...")

# Load latent and filtered metadata
latent = np.load(os.path.expanduser("~/Downloads/scvi_latent.npy"))
metadata = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata_filtered.csv"),
    index_col=0)

print(f"Latent shape: {latent.shape}")
print(f"Metadata shape: {metadata.shape}")
print(f"\nDisease groups:\n{metadata['disease_group'].value_counts()}")

# Compute UMAP — n_jobs=1 critical for Mac
print("\nComputing UMAP (this takes 2-3 min)...")
import umap
reducer = umap.UMAP(
    n_components=2,
    random_state=42,
    n_neighbors=15,
    min_dist=0.3,
    n_jobs=1        # CRITICAL — fixes Mac segfault
)
embedding = reducer.fit_transform(latent)
np.save(os.path.expanduser("~/Downloads/scvi_umap.npy"), embedding)
print(f"UMAP done! Shape: {embedding.shape}")

# Plot 1 — Disease group
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

colors = {
    "Reference":  "grey",
    "Control":    "steelblue",
    "AD_Low_Int": "orange",
    "AD_High":    "red"
}

for group, color in colors.items():
    mask = (metadata["disease_group"] == group).values
    if mask.sum() > 0:
        axes[0,0].scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=color, label=f"{group} (n={mask.sum()})",
            s=1, alpha=0.5)
axes[0,0].legend(markerscale=5, fontsize=9)
axes[0,0].set_title("Disease Group — scVI UMAP", fontsize=13)
axes[0,0].set_xlabel("UMAP 1")
axes[0,0].set_ylabel("UMAP 2")

# Plot 2 — Braak stage
braak_colors = {
    "Reference": "grey",   "Braak 0": "#2166ac",
    "Braak II":  "#74add1", "Braak III": "#fee090",
    "Braak IV":  "#fdae61", "Braak V":   "#f46d43",
    "Braak VI":  "#d73027"
}
for stage, color in braak_colors.items():
    mask = (metadata["Braak.stage"] == stage).values
    if mask.sum() > 0:
        axes[0,1].scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=color, label=stage, s=1, alpha=0.5)
axes[0,1].legend(markerscale=5, fontsize=8)
axes[0,1].set_title("Braak Stage — scVI UMAP", fontsize=13)
axes[0,1].set_xlabel("UMAP 1")

# Plot 3 — ADNC severity
adnc_colors = {
    "Reference":    "grey",    "Not AD":       "#2166ac",
    "Low":          "#74add1", "Intermediate": "#fdae61",
    "High":         "#d73027"
}
for adnc, color in adnc_colors.items():
    mask = (metadata["ADNC"] == adnc).values
    if mask.sum() > 0:
        axes[1,0].scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=color, label=adnc, s=1, alpha=0.5)
axes[1,0].legend(markerscale=5, fontsize=9)
axes[1,0].set_title("ADNC Severity — scVI UMAP", fontsize=13)
axes[1,0].set_xlabel("UMAP 1")
axes[1,0].set_ylabel("UMAP 2")

# Plot 4 — APOE4 status
apoe_colors = {"Reference": "grey", "N": "steelblue", "Y": "red"}
for apoe, color in apoe_colors.items():
    mask = (metadata["APOE4.status"] == apoe).values
    if mask.sum() > 0:
        axes[1,1].scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=color, label=f"APOE4={apoe} (n={mask.sum()})",
            s=1, alpha=0.5)
axes[1,1].legend(markerscale=5, fontsize=9)
axes[1,1].set_title("APOE4 Status — scVI UMAP", fontsize=13)
axes[1,1].set_xlabel("UMAP 1")

plt.suptitle("SEA-AD Human Microglia — scVI Integration\nRavi Shanker Pandey, Ph.D.",
             fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(
    os.path.expanduser("~/Desktop/sea_ad_analysis/results/scVI_UMAP.png"),
    dpi=300, bbox_inches="tight")
plt.close()
print("UMAP plot saved!")

# Summary statistics
print("\n=== Summary ===")
print(f"Total cells: {len(metadata)}")
print(f"\nDisease distribution:")
print(metadata["disease_group"].value_counts())
print(f"\nBraak distribution:")
print(metadata["Braak.stage"].value_counts())
print("\nDone! Results in ~/Desktop/sea_ad_analysis/results/")
