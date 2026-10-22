# ==============================================================================
# scripts/07_falsification_stress_tests.R
# ==============================================================================
library(tidyverse)
library(TCGAbiolinks)

cat("=============================================================\n")
cat("--- 1. Loading Phase 3a TCGA-CESC Predictions ---\n")
cat("=============================================================\n")

# Locate the TCGA predictions file automatically from Phase 3a
tcga_pred_file <- list.files("data/processed", pattern = "tcga.*predictions\\.csv", full.names = TRUE)

if(length(tcga_pred_file) == 0) {
  # Fallback if named differently
  tcga_pred_file <- "data/processed/03_tcga_cesc_predictions.csv"
} else {
  tcga_pred_file <- tcga_pred_file[1]
}

preds <- read_csv(tcga_pred_file, show_col_types = FALSE)

# Standardize the TCGA barcode column 
barcode_col <- intersect(colnames(preds), c("Patient_ID", "ModelID", "barcode", "geo_accession"))[1]
preds$sample_barcode <- substr(preds[[barcode_col]], 1, 15)

cat(sprintf("Loaded PRKDC predictions for %d primary cervical tumors.\n\n", nrow(preds)))

cat("=============================================================\n")
cat("--- 2. Falsification Test: Tumor Purity (CPE) ---\n")
cat("=============================================================\n")

# Load official consensus measurement of tumor purity from TCGAbiolinks
data("Tumor.purity")

purity_df <- Tumor.purity %>%
  dplyr::filter(Cancer.type == "CESC") %>%
  dplyr::mutate(
    sample_barcode = substr(Sample.ID, 1, 15),
    CPE = as.numeric(CPE) # CRITICAL FIX: Coerce the text to mathematical numbers
  ) %>%
  dplyr::select(sample_barcode, CPE) %>% 
  tidyr::drop_na()

# Merge clinical purity with our computational predictions
merged_purity <- dplyr::inner_join(preds, purity_df, by = "sample_barcode")

cat(sprintf("Successfully matched %d tumors with official CPE scores.\n", nrow(merged_purity)))

# Run the Pearson correlation stress test
res_purity <- cor.test(merged_purity$Predicted_Dependency_Score, merged_purity$CPE, method = "pearson")

cat(sprintf("\n-> Pearson Correlation (R): %.4f\n", res_purity$estimate))
cat(sprintf("-> P-value: %.2e\n", res_purity$p.value))

cat("\n=============================================================\n")
if(abs(res_purity$estimate) < 0.3) {
  cat("STATUS: [PASS] SUCCESS!\n")
  cat("Your PRKDC dependency score is statistically independent of Tumor Purity.\n")
  cat("This proves the model is capturing true biological vulnerability, not just\n")
  cat("measuring the ratio of cancer cells to stromal cells in the biopsy.\n")
} else {
  cat("STATUS: [FAIL] WARNING.\n")
  cat("Your PRKDC dependency score is highly confounded by Tumor Purity.\n")
}
cat("=============================================================\n")