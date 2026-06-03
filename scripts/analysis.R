# Save your complete analysis script

# SEA-AD Human Microglia Analysis
# Ravi Shanker Pandey, Ph.D. | The Jackson Laboratory
# Cross-species validation of AD microglial signatures

# Load libraries
library(Seurat)
library(zellkonverter)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(org.Hs.eg.db)

# Load SEA-AD data
sce <- readH5AD("data/microglia_dlpfc.h5ad")

# Check what's inside
cat("Dimensions:", dim(sce), "\n")
cat("\nColumn names (metadata):\n")
print(colnames(colData(sce)))

# Step 2 — Explore the metadata
# Check disease status
cat("=== Cognitive Status ===\n")
print(table(sce$Cognitive.status))

cat("\n=== ADNC severity ===\n")
print(table(sce$ADNC))

cat("\n=== Braak stages ===\n")
print(table(sce$Braak.stage))

cat("\n=== APOE4 status ===\n")
print(table(sce$APOE4.status))

cat("\n=== Microglia Subtypes ===\n")
print(table(sce$Subclass))

cat("\n=== Microglia Supertypes ===\n")
print(table(sce$Supertype))

# Convert to Seurat
seurat_mg <- as.Seurat(sce, counts = "X", data = NULL)
DefaultAssay(seurat_mg) <- "originalexp"

# Add disease groups
seurat_mg$disease_group <- case_when(
    seurat_mg$Cognitive.status == "Reference" ~ "Reference",
    seurat_mg$ADNC == "Not AD" ~ "Control",
    seurat_mg$ADNC %in% c("Low", "Intermediate") ~ "AD_Low_Int",
    seurat_mg$ADNC == "High" ~ "AD_High",
    TRUE ~ "Other"
)

# Normalize and cluster
seurat_mg <- NormalizeData(seurat_mg)
seurat_mg <- FindVariableFeatures(seurat_mg, nfeatures = 2000)
seurat_mg <- ScaleData(seurat_mg)
seurat_mg <- RunPCA(seurat_mg, npcs = 30)
seurat_mg <- RunUMAP(seurat_mg, dims = 1:15)
seurat_mg <- FindNeighbors(seurat_mg, dims = 1:15)
seurat_mg <- FindClusters(seurat_mg, resolution = 0.3)


# Plot 1 — UMAP colored by cluster
p1 <- DimPlot(seurat_mg, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              label = TRUE,
              repel = TRUE) +
  ggtitle("Microglia Clusters") +
  theme_minimal()

# Plot 2 — UMAP colored by disease status
p2 <- DimPlot(seurat_mg,
              reduction = "umap",
              group.by = "disease_group",
              cols = c("Reference" = "grey80",
                       "Control" = "steelblue",
                       "AD_Low_Int" = "orange",
                       "AD_High" = "red3")) +
  ggtitle("Disease Group") +
  theme_minimal()

# Plot 3 — UMAP colored by SEA-AD supertype
p3 <- DimPlot(seurat_mg,
              reduction = "umap",
              group.by = "Supertype",
              label = TRUE,
              repel = TRUE) +
  ggtitle("SEA-AD Supertypes") +
  theme_minimal()

# Plot 4 — UMAP colored by Braak stage
p4 <- DimPlot(seurat_mg,
              reduction = "umap",
              group.by = "Braak.stage") +
  ggtitle("Braak Stage") +
  theme_minimal()

# Combine plots
(p1 | p2) / (p3 | p4)

ggsave("plots/SEA_AD_microglia_UMAP.png",
       width = 16, height = 12, dpi = 300)


# Find DAM markers
# Focus on clusters 8, 9, 10, 12 — likely DAM states
dam_markers <- FindMarkers(seurat_mg,
                           ident.1 = c("8","9","10","12"),
                           ident.2 = c("0","1","2","3"),
                           min.pct = 0.25,
                           logfc.threshold = 0.5)

# Map to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(dam_markers),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add symbols to markers table
dam_markers$gene_symbol <- gene_symbols

# Show top markers with symbols
cat("=== Top DAM markers with gene symbols ===\n")
top_markers <- dam_markers[order(dam_markers$avg_log2FC, 
                                 decreasing=TRUE), ]
print(head(top_markers[, c("gene_symbol", "avg_log2FC", 
                           "pct.1", "pct.2", "p_val_adj")], 20))


# Convert all gene names to symbols
all_symbols <- mapIds(org.Hs.eg.db,
                      keys = rownames(seurat_mg),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

gene_map <- data.frame(ensembl = names(all_symbols),
                       symbol = all_symbols,
                       stringsAsFactors = FALSE)

# Key AD genes
key_genes <- c("TREM2","TYROBP","SPP1","P2RY12",
               "CX3CR1","CSF1R","TMEM119")
key_ensembl <- gene_map[gene_map$symbol %in% key_genes,]

# Create named vector for renaming
gene_name_map <- setNames(key_ensembl$symbol, 
                          key_ensembl$ensembl)

# Plot with gene symbols as titles
plots <- VlnPlot(seurat_mg,
                 features = key_ensembl$ensembl,
                 group.by = "disease_group",
                 ncol = 3,
                 cols = c("Reference" = "grey80",
                          "Control" = "steelblue",
                          "AD_Low_Int" = "orange",
                          "AD_High" = "red3"),
                 combine = FALSE)

# Rename each plot title
for(i in seq_along(plots)) {
  ensembl_id <- key_ensembl$ensembl[i]
  symbol <- key_ensembl$symbol[i]
  plots[[i]] <- plots[[i]] + 
    ggtitle(symbol) +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(face="bold", size=14))
}

# Combine and save
combined <- wrap_plots(plots, ncol=3)
ggsave("plots/SEA_AD_key_genes_symbols.png",
       plot = combined,
       width = 18, height = 12, dpi = 300)


# Expression trajectory
gene_summary <- data.frame()
for(i in 1:nrow(key_ensembl)) {
    ensembl <- key_ensembl$ensembl[i]
    symbol <- key_ensembl$symbol[i]
    expr <- FetchData(seurat_mg, vars = c(ensembl, "disease_group"))
    means <- expr %>%
        group_by(disease_group) %>%
        summarise(mean_expr = mean(.data[[ensembl]], na.rm=TRUE)) %>%
        mutate(gene = symbol)
    gene_summary <- rbind(gene_summary, means)
}


# Plot expression trend across disease severity
gene_summary$disease_group <- factor(
  gene_summary$disease_group,
  levels = c("Reference", "Control", "AD_Low_Int", "AD_High")
)

ggplot(gene_summary, 
       aes(x = disease_group, 
           y = mean_expr, 
           color = gene,
           group = gene)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title = "AD Gene Expression Across Disease Severity",
       subtitle = "SEA-AD Human Microglia DLPFC",
       x = "Disease Group",
       y = "Mean Expression",
       color = "Gene") +
  scale_color_brewer(palette = "Set1")

ggsave("plots/SEA_AD_disease_trajectory.png",
       width = 10, height = 6, dpi = 300)

cat("Trajectory plot saved!\n")



# Create publication-quality summary plot
# Separate into homeostatic vs DAM genes

gene_summary$gene_type <- case_when(
  gene_summary$gene %in% c("P2RY12", "CX3CR1", "TMEM119") ~ "Homeostatic",
  gene_summary$gene %in% c("SPP1", "TYROBP", "TREM2", "CSF1R") ~ "DAM/Activated",
  TRUE ~ "Other"
)

p_final <- ggplot(gene_summary,
                  aes(x = disease_group,
                      y = mean_expr,
                      color = gene,
                      group = gene)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~gene_type, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        strip.text = element_text(face="bold", size=12),
        legend.position = "right") +
  labs(title = "Homeostatic vs DAM Gene Expression in Human AD Microglia",
       subtitle = "SEA-AD DLPFC — Cross-species validation of mouse model findings",
       x = "Disease Group",
       y = "Mean Expression",
       color = "Gene") +
  scale_color_brewer(palette = "Set1")

ggsave("plots/SEA_AD_homeostatic_vs_DAM.png",
       plot = p_final,
       width = 14, height = 6, dpi = 300)

cat("Final plot saved!\n")


# Save results
saveRDS(seurat_mg, "results/SEA_AD_microglia_analyzed.rds")
write.csv(gene_summary, "results/gene_expression_summary.csv",
          row.names=FALSE)
write.csv(dam_markers, "results/DAM_markers.csv")
