
library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load full clustered metadata
meta_full <- read.csv(
  "~/Downloads/sea_ad_metadata_full_clustered.csv",
  row.names=1)

cat("Full metadata:", nrow(meta_full), "cells\n")
cat("Clusters:\n")
print(table(meta_full$scVI_clusters_full))

# Save full UMAP from Python
# Run in terminal first:
# python3 -c "
# import numpy as np, pandas as pd, os
# emb = np.load(os.path.expanduser('~/Downloads/scvi_umap_full.npy'))
# meta = pd.read_csv(os.path.expanduser('~/Downloads/sea_ad_metadata_full_clustered.csv'), index_col=0)
# df = pd.DataFrame(emb[:len(meta)], columns=['UMAP1','UMAP2'], index=meta.index)
# df.to_csv(os.path.expanduser('~/Downloads/scvi_umap_full_coords.csv'))
# print('Saved:', df.shape)
# "

# Load UMAP coordinates
umap_coords <- read.csv(
  "~/Downloads/scvi_umap_full_coords.csv",
  row.names=1)

cat("UMAP coords:", nrow(umap_coords), "cells\n")

# Match cells
common_cells <- intersect(rownames(meta_full),
                          rownames(umap_coords))
cat("Common cells:", length(common_cells), "\n")

meta_full <- meta_full[common_cells,]
umap_mat <- as.matrix(umap_coords[common_cells,])

# Create SCE
sce_full <- SingleCellExperiment(
  colData=meta_full)
reducedDim(sce_full, "UMAP") <- umap_mat
sce_full$scVI_clusters_full <- factor(
  meta_full$scVI_clusters_full)

cat("SCE ready:", ncol(sce_full), "cells\n")

# Run Slingshot
# Start = Cluster 4 (most homeostatic — highest Reference %)
# End = Cluster 0 (most DAM — 72% AD_High)
cat("Running Slingshot...\n")
sce_full <- slingshot(
  sce_full,
  clusterLabels = "scVI_clusters_full",
  reducedDim    = "UMAP",
  start.clus    = "4",
  end.clus      = "0"
)

cat("Slingshot complete!\n")
cat("Lineages:\n")
print(slingLineages(sce_full))

# Get pseudotime
pt_full <- slingPseudotime(sce_full)
cat("\nPseudotime summary:\n")
print(summary(pt_full))


# Get pseudotime
pt1 <- pt_full[,1]  # Lineage 1
pt2 <- pt_full[,2]  # Lineage 2

# Colors
colors_pt <- colorRampPalette(
  c("steelblue", "yellow", "red"))(100)

# ── Plot 1 — Trajectory UMAP ─────────────────────────────────────────────────
png("~/Desktop/sea_ad_analysis/plots/slingshot_full_trajectory.png",
    width=3000, height=900, res=150)
par(mfrow=c(1,3))

# Panel 1 — Lineage 1 pseudotime
pt1_scaled <- (pt1 - min(pt1, na.rm=TRUE)) /
  (max(pt1, na.rm=TRUE) - min(pt1, na.rm=TRUE))
pt1_scaled[is.na(pt1_scaled)] <- 0
col_idx <- pmax(1, pmin(100, cut(pt1_scaled,
                                 breaks=100, labels=FALSE)))
col_idx[is.na(col_idx)] <- 1

plot(umap_mat[,1], umap_mat[,2],
     col=colors_pt[col_idx],
     pch=16, cex=0.2,
     main="Lineage 1 Pseudotime\n4→2→5→0 (Primary DAM)",
     xlab="UMAP 1", ylab="UMAP 2")
lines(slingCurves(sce_full)[[1]]$s[
  slingCurves(sce_full)[[1]]$ord,],
  lwd=3, col="black")

# Panel 2 — Lineage 2 pseudotime
pt2_scaled <- (pt2 - min(pt2, na.rm=TRUE)) /
  (max(pt2, na.rm=TRUE) - min(pt2, na.rm=TRUE))
pt2_scaled[is.na(pt2_scaled)] <- 0
col_idx2 <- pmax(1, pmin(100, cut(pt2_scaled,
                                  breaks=100, labels=FALSE)))
col_idx2[is.na(col_idx2)] <- 1

plot(umap_mat[,1], umap_mat[,2],
     col=colors_pt[col_idx2],
     pch=16, cex=0.2,
     main="Lineage 2 Pseudotime\n4→2→3→1 (Alternative)",
     xlab="UMAP 1", ylab="UMAP 2")
lines(slingCurves(sce_full)[[2]]$s[
  slingCurves(sce_full)[[2]]$ord,],
  lwd=3, col="black")

# Panel 3 — Disease group
disease_cols <- c(
  "Reference"="grey",  "Control"="steelblue",
  "AD_Low_Int"="orange", "AD_High"="red")
plot(umap_mat[,1], umap_mat[,2],
     col=disease_cols[meta_full$disease_group],
     pch=16, cex=0.2,
     main="Disease Group",
     xlab="UMAP 1", ylab="UMAP 2")
legend("topright", legend=names(disease_cols),
       col=disease_cols, pch=16, cex=0.8)
lines(slingCurves(sce_full)[[1]]$s[
  slingCurves(sce_full)[[1]]$ord,],
  lwd=3, col="black")
lines(slingCurves(sce_full)[[2]]$s[
  slingCurves(sce_full)[[2]]$ord,],
  lwd=3, col="darkgreen")

dev.off()
cat("Trajectory plot saved!\n")

# ── Plot 2 — Pseudotime violin ───────────────────────────────────────────────
pt_df <- data.frame(
  pseudotime_L1 = pt1,
  pseudotime_L2 = pt2,
  disease_group = meta_full$disease_group,
  cluster       = as.character(
    meta_full$scVI_clusters_full),
  Braak         = meta_full$Braak.stage,
  ADNC          = meta_full$ADNC
) %>% tidyr::drop_na(pseudotime_L1)

p1 <- ggplot(pt_df,
             aes(x=factor(disease_group,
                          levels=c("Reference","Control",
                                   "AD_Low_Int","AD_High")),
                 y=pseudotime_L1,
                 fill=disease_group)) +
  geom_violin(alpha=0.7) +
  geom_boxplot(width=0.1, fill="white",
               outlier.size=0.1) +
  scale_fill_manual(values=c(
    "Reference"="grey",   "Control"="steelblue",
    "AD_Low_Int"="orange","AD_High"="red")) +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Lineage 1 Pseudotime by Disease Group",
       subtitle="Primary DAM path: Cluster 4→2→5→0",
       x="Disease Group", y="Pseudotime")

# Braak stage violin
pt_df_braak <- pt_df %>%
  filter(!is.na(pseudotime_L1)) %>%
  filter(Braak %in% c("Braak 0","Braak II",
                      "Braak III","Braak IV",
                      "Braak V","Braak VI"))

p2 <- ggplot(pt_df_braak,
             aes(x=factor(Braak,
                          levels=c("Braak 0","Braak II",
                                   "Braak III","Braak IV",
                                   "Braak V","Braak VI")),
                 y=pseudotime_L1,
                 fill=Braak)) +
  geom_violin(alpha=0.7) +
  geom_boxplot(width=0.1, fill="white",
               outlier.size=0.1) +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Lineage 1 Pseudotime by Braak Stage",
       subtitle="Progressive increase with tau pathology",
       x="Braak Stage", y="Pseudotime")

library(patchwork)
combined <- p1 + p2
ggsave("~/Desktop/sea_ad_analysis/plots/slingshot_full_violin.png",
       plot=combined, width=14, height=6, dpi=300)
cat("Violin plots saved!\n")

# ── Summary statistics ────────────────────────────────────────────────────────
cat("\n=== Mean Lineage 1 pseudotime by disease ===\n")
print(tapply(pt1, meta_full$disease_group,
             mean, na.rm=TRUE))

cat("\n=== Mean Lineage 1 pseudotime by cluster ===\n")
print(tapply(pt1, meta_full$scVI_clusters_full,
             mean, na.rm=TRUE))

cat("\n=== Mean Lineage 1 pseudotime by Braak ===\n")
print(tapply(pt1, meta_full$Braak.stage,
             mean, na.rm=TRUE))

# Save results
write.csv(pt_df,
          "~/Desktop/sea_ad_analysis/results/slingshot_full_pseudotime.csv",
          row.names=FALSE)
cat("\nResults saved!\n")