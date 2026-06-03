# Add to scvi_markers_fix.py after markers are extracted
# Save as ~/sea_ad_analysis/scripts/convert_markers.py

import pandas as pd
import os

# Load org.Hs.eg.db mapping from R
# First run this in R:
# library(org.Hs.eg.db)
# gene_map <- select(org.Hs.eg.db,
#                    keys=keys(org.Hs.eg.db, "ENSEMBL"),
#                    columns=c("ENSEMBL","SYMBOL"),
#                    keytype="ENSEMBL")
# write.csv(gene_map, "~/Downloads/human_gene_map.csv", row.names=FALSE)

gene_map = pd.read_csv(
    os.path.expanduser("~/Downloads/human_gene_map.csv"))
gene_map = gene_map.dropna().drop_duplicates("ENSEMBL")
ensembl_to_symbol = dict(zip(gene_map["ENSEMBL"], gene_map["SYMBOL"]))

print(f"Gene map loaded: {len(ensembl_to_symbol)} genes")

# Convert markers for each cluster
for cluster in ["0", "1", "2", "3", "4"]:
    df = pd.read_csv(os.path.expanduser(
        f"~/Desktop/sea_ad_analysis/results/markers/cluster_{cluster}_markers.csv"))
    
    df["symbol"] = df["gene"].map(ensembl_to_symbol)
    df = df[["symbol", "gene", "logfoldchange", "pval_adj", "score"]]
    df = df.sort_values("logfoldchange", ascending=False)
    
    print(f"\nCluster {cluster} top 10 (with symbols):")
    print(df[["symbol","logfoldchange"]].head(10).to_string())
    
    df.to_csv(os.path.expanduser(
        f"~/Desktop/sea_ad_analysis/results/markers/cluster_{cluster}_markers_annotated.csv"),
        index=False)

print("\nAnnotated markers saved!")
