# Save as ~/Desktop/sea_ad_analysis/scripts/scvi_umap_full.py
# Then run:
import os
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap

print("Loading full dataset results...")

latent = np.load(os.path.expanduser("~/Downloads/scvi_latent_full.npy"))
metadata = pd.read_csv(
    os.path.expanduser("~/Downloads/sea_ad_metadata_full_filtered.csv"),
    index_col=0)

print(f"Latent: {latent.shape}")
print(f"Metadata: {metadata.shape}")

# Compute UMAP
print("\nComputing UMAP (5-10 min)...")
reducer = umap.UMAP(
    n_components=2,
    random_state=42,
    n_neighbors=15,
    min_dist=0.3,
    n_jobs=1
)
embedding = reducer.fit_transform(latent)
np.save(os.path.expanduser("~/Downloads/scvi_umap_full.npy"), embedding)
print(f"UMAP done! Shape: {embedding.shape}")

# Plot
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

disease_colors = {
    "Reference":  "grey",
    "Control":    "steelblue",
    "AD_Low_Int": "orange",
    "AD_High":    "red"
}

for group, color in disease_colors.items():
    mask = (metadata["disease_group"] == group).values
    axes[0,0].scatter(embedding[mask,0], embedding[mask,1],
                     c=color, label=f"{group} (n={mask.sum()})",
                     s=0.5, alpha=0.4)
axes[0,0].legend(markerscale=6, fontsize=9)
axes[0,0].set_title("Disease Group — Full Dataset scVI", fontsize=13, fontweight="bold")
axes[0,0].set_xlabel("UMAP 1")
axes[0,0].set_ylabel("UMAP 2")

braak_colors = {
    "Reference": "grey",   "Braak 0":  "#2166ac",
    "Braak II":  "#74add1", "Braak III":"#fee090",
    "Braak IV":  "#fdae61", "Braak V":  "#f46d43",
    "Braak VI":  "#d73027"
}
for stage, color in braak_colors.items():
    mask = (metadata["Braak.stage"] == stage).values
    if mask.sum() > 0:
        axes[0,1].scatter(embedding[mask,0], embedding[mask,1],
                         c=color, label=stage, s=0.5, alpha=0.4)
axes[0,1].legend(markerscale=6, fontsize=8)
axes[0,1].set_title("Braak Stage", fontsize=13, fontweight="bold")
axes[0,1].set_xlabel("UMAP 1")

adnc_colors = {
    "Reference": "grey",    "Not AD":       "#2166ac",
    "Low":       "#74add1", "Intermediate": "#fdae61",
    "High":      "#d73027"
}
for adnc, color in adnc_colors.items():
    mask = (metadata["ADNC"] == adnc).values
    if mask.sum() > 0:
        axes[1,0].scatter(embedding[mask,0], embedding[mask,1],
                         c=color, label=adnc, s=0.5, alpha=0.4)
axes[1,0].legend(markerscale=6, fontsize=9)
axes[1,0].set_title("ADNC Severity", fontsize=13, fontweight="bold")
axes[1,0].set_xlabel("UMAP 1")
axes[1,0].set_ylabel("UMAP 2")

apoe_colors = {"Reference": "grey", "N": "steelblue", "Y": "red"}
for apoe, color in apoe_colors.items():
    mask = (metadata["APOE4.status"] == apoe).values
    if mask.sum() > 0:
        axes[1,1].scatter(embedding[mask,0], embedding[mask,1],
                         c=color,
                         label=f"APOE4={apoe} (n={mask.sum()})",
                         s=0.5, alpha=0.4)
axes[1,1].legend(markerscale=6, fontsize=9)
axes[1,1].set_title("APOE4 Status", fontsize=13, fontweight="bold")
axes[1,1].set_xlabel("UMAP 1")

plt.suptitle("SEA-AD Human Microglia — Full Dataset scVI\n(42,486 cells, 42 donors)",
             fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(os.path.expanduser(
    "~/Desktop/sea_ad_analysis/results/scVI_UMAP_full.png"),
    dpi=300, bbox_inches="tight")
plt.close()
print("UMAP plot saved!")
