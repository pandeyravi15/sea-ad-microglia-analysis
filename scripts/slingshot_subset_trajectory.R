library(slingshot)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)

# Load your existing Seurat microglia object
seurat_mg <- readRDS("~/Downloads/SEA_AD_microglia_analyzed.rds")

# Add scVI clusters to Seurat object
scvi_meta <- read.csv(
  "~/Downloads/sea_ad_metadata_clustered.csv",
  row.names = 1)

# Add scVI clusters
seurat_mg$scVI_clusters <- scvi_meta[
  colnames(seurat_mg), "scVI_clusters"]

cat("Seurat object loaded!\n")
cat("Cells:", ncol(seurat_mg), "\n")
cat("scVI clusters:\n")
print(table(seurat_mg$scVI_clusters))



# ── Step 1: Convert Seurat to SingleCellExperiment ───────────────────────────
cat("Converting to SingleCellExperiment...\n")

sce <- as.SingleCellExperiment(seurat_mg)

# Add scVI UMAP coordinates from Python results
umap_coords <- read.csv(
  "~/Downloads/sea_ad_metadata_clustered.csv",
  row.names = 1)[, c("scVI_clusters")]

# Load scVI UMAP embedding
# First save from Python — run this in terminal:
# python3 -c "
# import numpy as np, pandas as pd
# emb = np.load('/Users/pandera/Downloads/scvi_umap.npy')
# meta = pd.read_csv('/Users/pandera/Downloads/sea_ad_metadata_clustered.csv', index_col=0)
# df = pd.DataFrame(emb[:len(meta)], columns=['UMAP1','UMAP2'], index=meta.index)
# df.to_csv('/Users/pandera/Downloads/scvi_umap_coords.csv')
# print('Saved!')
# "

umap_df <- read.csv(
  "~/Downloads/scvi_umap_coords.csv",
  row.names = 1)

cat("UMAP coords loaded:", nrow(umap_df), "cells\n")

# Match cells
common_cells <- intersect(colnames(sce), rownames(umap_df))
cat("Common cells:", length(common_cells), "\n")

sce <- sce[, common_cells]
umap_matrix <- as.matrix(umap_df[common_cells, ])
reducedDim(sce, "UMAP") <- umap_matrix

# Add scVI clusters
sce$scVI_clusters <- factor(
  scvi_meta[common_cells, "scVI_clusters"])

cat("SCE ready:", ncol(sce), "cells\n")


# ── Step 2: Run Slingshot ────────────────────────────────────────────────────
cat("Running Slingshot trajectory analysis...\n")

# Run slingshot
# Start cluster = Cluster 2 (most homeostatic — highest P2RY12/CX3CR1)
# End cluster = Cluster 1 (most DAM — 67% AD_High)
sce <- slingshot(
  sce,
  clusterLabels = "scVI_clusters",
  reducedDim    = "UMAP",
  start.clus    = "2",    # homeostatic start
  end.clus      = "1"     # DAM end
)

cat("Slingshot complete!\n")
cat("Lineages found:\n")
print(slingLineages(sce))

# Get pseudotime
pseudotime <- slingPseudotime(sce)
cat("\nPseudotime summary:\n")
print(summary(pseudotime))


# ── Step 3: Visualize trajectories ───────────────────────────────────────────

# Get UMAP coordinates
umap_plot <- reducedDim(sce, "UMAP")

# Get pseudotime for lineage 1
pt <- pseudotime[, 1]
pt_scaled <- (pt - min(pt, na.rm=TRUE)) /
  (max(pt, na.rm=TRUE) - min(pt, na.rm=TRUE))

# Color palette
colors_pt <- colorRampPalette(
  c("steelblue", "yellow", "red"))(100)

# ── Plot 1 — Pseudotime on UMAP ──────────────────────────────────────────────
png("~/Desktop/sea_ad_analysis/plots/slingshot_pseudotime.png",
    width=2400, height=800, res=150)

par(mfrow=c(1,3))

# Panel 1 — Pseudotime
color_idx <- cut(pt_scaled, breaks=100, labels=FALSE)
color_idx[is.na(color_idx)] <- 1

plot(umap_plot[,1], umap_plot[,2],
     col  = colors_pt[color_idx],
     pch  = 16, cex = 0.3,
     main = "Pseudotime\n(Blue=early → Red=late)",
     xlab = "UMAP 1", ylab = "UMAP 2")

# Add trajectory curves
lines(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord, ],
      lwd=3, col="black")

# Panel 2 — Clusters
cluster_colors <- c("0"="steelblue", "1"="red",
                    "2"="green3",    "3"="orange",
                    "4"="purple")
plot(umap_plot[,1], umap_plot[,2],
     col  = cluster_colors[as.character(sce$scVI_clusters)],
     pch  = 16, cex = 0.3,
     main = "scVI Clusters",
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("topright",
       legend = names(cluster_colors),
       col    = cluster_colors,
       pch    = 16, cex = 0.8,
       title  = "Cluster")
lines(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord,],
      lwd=3, col="black")

# Panel 3 — Disease group
disease_cols <- c(
  "Reference"  = "grey",
  "Control"    = "steelblue",
  "AD_Low_Int" = "orange",
  "AD_High"    = "red"
)
plot(umap_plot[,1], umap_plot[,2],
     col  = disease_cols[sce$disease_group],
     pch  = 16, cex = 0.3,
     main = "Disease Group",
     xlab = "UMAP 1", ylab = "UMAP 2")
legend("topright",
       legend = names(disease_cols),
       col    = disease_cols,
       pch    = 16, cex = 0.8)
lines(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord,],
      lwd=3, col="black")

dev.off()
cat("Trajectory plot saved!\n")

# ── Plot 2 — Pseudotime by disease group ─────────────────────────────────────
pt_df <- data.frame(
  pseudotime    = pt,
  disease_group = sce$disease_group,
  cluster       = as.character(sce$scVI_clusters)
) %>% tidyr::drop_na()

p <- ggplot(pt_df, aes(x = factor(disease_group,
                                  levels = c("Reference","Control",
                                             "AD_Low_Int","AD_High")),
                       y = pseudotime,
                       fill = disease_group)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white",
               outlier.size = 0.1) +
  scale_fill_manual(values = c(
    "Reference"  = "grey",
    "Control"    = "steelblue",
    "AD_Low_Int" = "orange",
    "AD_High"    = "red")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title    = "Pseudotime by Disease Group",
       subtitle = "Slingshot trajectory — Homeostatic → DAM",
       x = "Disease Group",
       y = "Pseudotime")

ggsave("~/Desktop/sea_ad_analysis/plots/slingshot_pseudotime_violin.png",
       plot=p, width=8, height=6, dpi=300)
cat("Violin plot saved!\n")

# ── Summary statistics ────────────────────────────────────────────────────────
cat("\n=== Mean pseudotime by disease group ===\n")
print(tapply(pt, sce$disease_group, mean, na.rm=TRUE))

cat("\n=== Mean pseudotime by cluster ===\n")
print(tapply(pt, sce$scVI_clusters, mean, na.rm=TRUE))


# Save pseudotime results
pt_df <- data.frame(
  cell = colnames(sce),
  pseudotime_lineage1 = slingPseudotime(sce)[,1],
  pseudotime_lineage2 = slingPseudotime(sce)[,2],
  disease_group = sce$disease_group,
  scVI_cluster = sce$scVI_clusters,
  Braak_stage = sce$Braak.stage,
  ADNC = sce$ADNC
)

write.csv(pt_df,
          "~/Desktop/sea_ad_analysis/results/slingshot_pseudotime.csv",
          row.names = FALSE)

cat("Pseudotime results saved!\n")
cat("\nMean pseudotime by disease group:\n")
print(tapply(pt_df$pseudotime_lineage2,
             pt_df$disease_group, mean, na.rm=TRUE))
