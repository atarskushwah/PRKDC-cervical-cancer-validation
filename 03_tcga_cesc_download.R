# ==============================================================================
# scripts/03_tcga_cesc_projection.R
# ==============================================================================
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(tidymodels)

cat("--- 1. Data & Model Initialization ---\n")
locked_recipe <- readRDS("models/locked_preprocessing_recipe.rds")
locked_model <- readRDS("models/locked_prkdc_dependency_model.rds")

expected_vars <- locked_recipe$var_info$variable
template_data <- locked_recipe$template

cat("--- 2. Querying TCGA-CESC Data ---\n")
# Query standard STAR RNA-seq counts for Cervical Cancer
query <- GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download and parse data (This may take a moment depending on the GDC server)
GDCdownload(query)
tcga_data <- GDCprepare(query)

cat("--- 3. Processing Expression Data ---\n")
# Extract TPM (closest to DepMap transcriptomics scale) and clinical data
expr_tpm <- assay(tcga_data, "tpm_unstrand")
clinical <- as.data.frame(colData(tcga_data))

# Force Log2(TPM + 1) transform to identically match the DepMap training scale
cat("Applying log2(TPM + 1) transform...\n")
expr_log2 <- log2(expr_tpm + 1)

# Map Ensembl IDs to standard Gene Symbols
gene_info <- as.data.frame(rowData(tcga_data))
symbols <- gene_info$gene_name

valid_idx <- which(!is.na(symbols) & symbols != "")
expr_mapped <- expr_log2[valid_idx, , drop = FALSE]
mapped_symbols <- symbols[valid_idx]

cat("Aggregating duplicate genes (with NA protection)...\n")
expr_df <- as.data.frame(expr_mapped)
expr_df$Symbol <- mapped_symbols

expr_agg <- expr_df %>%
  group_by(Symbol) %>%
  summarise(across(everything(), ~ suppressWarnings(max(.x, na.rm = TRUE)))) %>%
  ungroup()

# Clean up any potential -Inf values and transpose for prediction
expr_agg_mat <- as.matrix(expr_agg[, -1])
expr_agg_mat[is.infinite(expr_agg_mat)] <- 0

tcga_expr_t <- as.data.frame(t(expr_agg_mat))
colnames(tcga_expr_t) <- expr_agg$Symbol
tcga_expr_t$barcode <- rownames(tcga_expr_t)

tcga_cohort <- tcga_expr_t
tcga_cohort$OncotreeLineage <- "Cervix"
tcga_cohort$ModelID <- tcga_cohort$barcode

cat("--- 4. Direct Type-Safe Data Injection ---\n")
current_cols <- colnames(tcga_cohort)
clean_current_cols <- gsub("[^A-Z0-9]", "", toupper(current_cols))
matched_count <- 0

for (expected_col in expected_vars) {
  if(expected_col %in% c("ModelID", "OncotreeLineage", "PrimaryOrMetastasis", "PRKDC_effect")) next
  if (endsWith(tolower(expected_col), "_cn") || endsWith(tolower(expected_col), "_mut")) next
  
  base_expected <- sub(" \\([0-9]+\\).*|\\.\\.[0-9]+\\..*", "", expected_col)
  clean_expected <- gsub("[^A-Z0-9]", "", toupper(base_expected))
  match_idx <- which(clean_current_cols == clean_expected)
  
  if (length(match_idx) > 0) {
    colnames(tcga_cohort)[match_idx[1]] <- expected_col
    clean_current_cols[match_idx[1]] <- "MAPPED" 
    matched_count <- matched_count + 1
  }
}
cat(sprintf("Successfully mapped %d biological features!\n", matched_count))

cat("Constructing perfectly typed dataframe...\n")
n_patients <- nrow(tcga_cohort)
aligned_list <- list()

for (col in expected_vars) {
  if (col %in% colnames(tcga_cohort)) {
    val <- tcga_cohort[[col]]
    if(is.data.frame(val) || is.matrix(val)) val <- val[[1]]
  } else {
    val <- rep(0, n_patients)
  }
  
  if (!is.null(template_data) && col %in% colnames(template_data)) {
    tmpl_class <- class(template_data[[col]])[1]
    
    if (tmpl_class == "numeric") val <- as.numeric(val)
    else if (tmpl_class == "integer") val <- as.integer(val)
    else if (tmpl_class == "logical") val <- as.logical(val)
    else if (tmpl_class == "character") val <- as.character(val)
    else if (tmpl_class == "factor") {
      val <- factor(as.character(val), levels = levels(template_data[[col]]))
    }
    
    if (tmpl_class %in% c("numeric", "integer")) val[is.na(val)] <- 0
    else if (tmpl_class == "logical") val[is.na(val)] <- FALSE
    else if (tmpl_class == "character") val[is.na(val)] <- "Unknown"
    else if (tmpl_class == "factor") val[is.na(val)] <- levels(template_data[[col]])[1]
    
  } else {
    val <- as.numeric(val)
    val[is.na(val)] <- 0
  }
  aligned_list[[col]] <- val
}

tcga_final <- as_tibble(aligned_list)
tcga_final$barcode <- tcga_cohort$barcode

cat("--- 5. Prediction & Repair State Mapping ---\n")
var_check <- apply(tcga_final %>% select_if(is.numeric), 2, var, na.rm = TRUE)
active_features <- sum(var_check > 0, na.rm = TRUE)
cat(sprintf("Matrix check: %d features have active biological variance.\n", active_features))

# Generate predictions for all TCGA tumors
preds <- predict(locked_model, new_data = tcga_final)

# Integrate predictions with TCGA's native clinical and HPV metadata
tcga_results <- tcga_final %>%
  select(barcode) %>%
  mutate(Predicted_Dependency_Score = preds[[1]]) %>%
  left_join(
    clinical %>% select(barcode, any_of(c("paper_HPV_Status", "primary_diagnosis", "tumor_stage"))), 
    by = "barcode"
  )

pred_sd <- sd(tcga_results$Predicted_Dependency_Score, na.rm = TRUE)
cat(sprintf("Prediction Variance (SD): %.5f\n", pred_sd))

write_csv(tcga_results, "data/processed/03_tcga_cesc_predictions.csv")
cat("\n=============================================================\n")
cat("SUCCESS: Saved TCGA-CESC predictions and clinical metadata to:\n")
cat("'data/processed/03_tcga_cesc_predictions.csv'\n")
cat("=============================================================\n")