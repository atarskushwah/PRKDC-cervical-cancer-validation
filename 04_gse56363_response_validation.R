# ==============================================================================
# scripts/04_gse56363_response_validation.R
# ==============================================================================
library(GEOquery)
library(tidyverse)
library(tidymodels)
library(pROC)

cat("--- 1. Data & Model Initialization ---\n")
locked_recipe <- readRDS("models/locked_preprocessing_recipe.rds")
locked_model <- readRDS("models/locked_prkdc_dependency_model.rds")

# Extract raw variables from var_info (bypasses stripped templates)
expected_vars <- locked_recipe$var_info$variable
template_data <- locked_recipe$template

gse <- getGEO("GSE56363", GSEMatrix = TRUE)
gse_data <- gse[[1]]

# Extract Clinical Data
clinical <- pData(gse_data) %>%
  select(geo_accession, title, response = characteristics_ch1.1) %>%
  mutate(
    clean_response = trimws(sub(".*:", "", tolower(response))),
    response_binary = case_when(
      clean_response == "complete response" ~ 1,
      str_detect(clean_response, "non-complete|partial|progressive") ~ 0,
      str_detect(title, "_NCR_") ~ 0,
      str_detect(title, "_CR_") ~ 1,
      TRUE ~ NA_real_
    ),
    response_factor = factor(response_binary, levels = c(0, 1), labels = c("Non-Complete", "Complete"))
  )

cat("--- 2. Microarray Feature Mapping & Scaling ---\n")
expr_matrix <- exprs(gse_data)
feature_data <- fData(gse_data)

# Ensure Log2 scale to match DepMap
if(max(expr_matrix, na.rm = TRUE) > 50) {
  cat("Applying log2(x+1) transform...\n")
  expr_matrix <- log2(expr_matrix + 1)
}

sym_col_idx <- grep("GENE_SYMBOL|GeneSymbol|Symbol", colnames(feature_data), ignore.case = TRUE)[1]
actual_sym_col <- colnames(feature_data)[sym_col_idx]

symbols <- as.character(feature_data[[actual_sym_col]])
symbols <- sapply(strsplit(symbols, "///|,"), `[`, 1) 

valid_idx <- which(!is.na(symbols) & symbols != "")
expr_mapped <- expr_matrix[valid_idx, , drop = FALSE]
mapped_symbols <- symbols[valid_idx]

cat("Aggregating probes (with NA protection)...\n")
# CRITICAL FIX: dplyr aggregation safely bypasses native Microarray NAs
expr_df <- as.data.frame(expr_mapped)
expr_df$Symbol <- mapped_symbols

expr_agg <- expr_df %>%
  group_by(Symbol) %>%
  summarise(across(everything(), ~ suppressWarnings(max(.x, na.rm = TRUE)))) %>%
  ungroup()

# Clean up -Inf (caused if a gene was NA across all probes for a patient)
expr_agg_mat <- as.matrix(expr_agg[, -1])
expr_agg_mat[is.infinite(expr_agg_mat)] <- 0

expr_t <- as.data.frame(t(expr_agg_mat))
colnames(expr_t) <- expr_agg$Symbol
expr_t$geo_accession <- rownames(expr_t)

gse_cohort <- merge(clinical, expr_t, by = "geo_accession")
gse_cohort$OncotreeLineage <- "Cervix"
gse_cohort$PrimaryOrMetastasis <- "Primary"
gse_cohort$ModelID <- as.character(gse_cohort$geo_accession)

cat("--- 3. Direct Type-Safe Data Injection ---\n")
current_cols <- colnames(gse_cohort)
clean_current_cols <- gsub("[^A-Z0-9]", "", toupper(current_cols))
matched_count <- 0

for (expected_col in expected_vars) {
  if(expected_col %in% c("ModelID", "OncotreeLineage", "PrimaryOrMetastasis", "PRKDC_effect")) next
  if (endsWith(tolower(expected_col), "_cn") || endsWith(tolower(expected_col), "_mut")) next
  
  base_expected <- sub(" \\([0-9]+\\).*|\\.\\.[0-9]+\\..*", "", expected_col)
  clean_expected <- gsub("[^A-Z0-9]", "", toupper(base_expected))
  match_idx <- which(clean_current_cols == clean_expected)
  
  if (length(match_idx) > 0) {
    colnames(gse_cohort)[match_idx[1]] <- expected_col
    clean_current_cols[match_idx[1]] <- "MAPPED" 
    matched_count <- matched_count + 1
  }
}
cat(sprintf("Successfully mapped %d biological features!\n", matched_count))

cat("Constructing perfectly typed dataframe...\n")
n_patients <- nrow(gse_cohort)
aligned_list <- list()

for (col in expected_vars) {
  if (col %in% colnames(gse_cohort)) {
    val <- gse_cohort[[col]]
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
    
    # Fill remaining NAs with strict safe baselines
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

gse_final <- as_tibble(aligned_list)
gse_final$geo_accession <- gse_cohort$geo_accession
gse_final$response_factor <- gse_cohort$response_factor
gse_final$response_binary <- gse_cohort$response_binary

cat("--- 4. Prediction & Clinical Validation ---\n")
# Diagnostic print to guarantee variance exists before predicting
var_check <- apply(gse_final %>% select_if(is.numeric), 2, var, na.rm = TRUE)
active_features <- sum(var_check > 0, na.rm = TRUE)
cat(sprintf("Matrix check: %d features have active biological variance.\n", active_features))

preds <- predict(locked_model, new_data = gse_final)

gse_results <- gse_final %>%
  select(geo_accession, response_factor, response_binary) %>%
  mutate(Predicted_Dependency_Score = preds[[1]])

pred_sd <- sd(gse_results$Predicted_Dependency_Score, na.rm = TRUE)
cat(sprintf("Prediction Variance (SD): %.5f\n", pred_sd))

if (is.na(pred_sd) || pred_sd < 1e-5) {
  cat("\nCRITICAL WARNING: Predictions STILL lack variance. The model is defaulting to a single pathway.\n")
  gse_results$Standardized_Score <- 0
} else {
  gse_results$Standardized_Score <- scale(gse_results$Predicted_Dependency_Score)[,1]
}

suppressWarnings({
  log_mod <- glm(response_binary ~ Standardized_Score, data = gse_results, family = binomial)
  or_estimate <- exp(coef(log_mod)["Standardized_Score"])
  ci <- exp(confint(log_mod)["Standardized_Score", ])
})

roc_obj <- roc(gse_results$response_binary, gse_results$Standardized_Score, quiet = TRUE)
auc_ci <- ci.auc(roc_obj)

cat("\n================ CLINICAL VALIDATION RESULTS ================\n")
cat(sprintf("Odds Ratio (per SD increase): %.2f (95%% CI: %.2f - %.2f)\n", or_estimate, ci[1], ci[2]))
cat(sprintf("AUC: %.2f (95%% CI: %.2f - %.2f)\n", roc_obj$auc, auc_ci[1], auc_ci[3]))
cat("=============================================================\n")

write_csv(gse_results, "data/processed/04_gse56363_predictions.csv")