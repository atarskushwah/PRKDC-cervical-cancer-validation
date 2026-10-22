# ==============================================================================
# scripts/12a_generate_pathways.R
# ==============================================================================
library(GEOquery)
library(dplyr)
library(tidyr)

setwd("~/CERV-DPK_Multiomics")

cat("Downloading GSE3578 and standardized annotations from GEO...\n")
# FIX: Added AnnotGPL = TRUE to guarantee standard Gene Symbol mappings
gse <- getGEO("GSE3578", GSEMatrix = TRUE, AnnotGPL = TRUE)
eset <- gse[[1]]

# Extract expression matrix and feature data
expr <- exprs(eset)
fdata <- fData(eset)

# Explicitly pull the standardized 'Gene symbol' column provided by AnnotGPL
if("Gene symbol" %in% colnames(fdata)) {
  # Split cases where multiple genes are mapped (e.g., "BRCA1///other")
  fdata$GeneSymbol <- sapply(strsplit(as.character(fdata$`Gene symbol`), "///"), `[`, 1)
} else {
  stop("Standardized Gene Symbol column still not found. Check network connection.")
}

# Average expression by Gene Symbol to handle multiple probes per gene
expr_df <- as.data.frame(expr)
expr_df$GeneSymbol <- trimws(fdata$GeneSymbol)

gene_expr <- expr_df %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

# Convert to matrix and calculate Z-scores per gene across the cohort
genes <- gene_expr$GeneSymbol
mat <- as.matrix(gene_expr[, -1])
rownames(mat) <- genes
mat_z <- t(scale(t(mat)))

# Define standard DNA repair signatures
hr_genes <- c("BRCA1", "BRCA2", "RAD51", "MRE11A", "NBN", "RAD50", "ATM")
nhej_genes <- c("PRKDC", "XRCC6", "XRCC5", "LIG4", "XRCC4")
rep_stress <- c("RPA1", "RPA2", "RPA3", "CHEK1")

# Function to safely calculate mean signature score for available genes
get_score <- function(sig) {
  avail <- intersect(sig, rownames(mat_z))
  if(length(avail) == 0) return(rep(NA, ncol(mat_z)))
  colMeans(mat_z[avail, , drop = FALSE], na.rm = TRUE)
}

# Calculate scores for each patient (geo_accession)
scores <- data.frame(
  geo_accession = colnames(mat_z),
  HR_Score = get_score(hr_genes),
  NHEJ_Score = get_score(nhej_genes),
  RepStress_Score = get_score(rep_stress)
)

# Merge with your existing predictions
pred <- read.csv("data/processed/05_gse3578_predictions.csv")
final_df <- merge(pred, scores, by = "geo_accession")

write.csv(final_df, "data/processed/05_gse3578_pathway_scores.csv", row.names = FALSE)
cat("SUCCESS: Pathway scores saved to 'data/processed/05_gse3578_pathway_scores.csv'!\n")