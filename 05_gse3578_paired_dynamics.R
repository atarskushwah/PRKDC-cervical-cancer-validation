# ==============================================================================
# scripts/05_gse3578_paired_dynamics.R
# ==============================================================================
library(GEOquery)
library(tidyverse)
library(tidymodels)
library(org.Hs.eg.db)

cat("--- 1. Data & Model Initialization ---\n")
locked_recipe <- readRDS("models/locked_preprocessing_recipe.rds")
locked_model <- readRDS("models/locked_prkdc_dependency_model.rds")

expected_vars <- locked_recipe$var_info$variable
template_data <- locked_recipe$template

cat("--- 2. Querying GSE3578 Paired Dynamics Cohort ---\n")
gse <- getGEO("GSE3578", GSEMatrix = TRUE)
gse_data <- gse[[1]]

clinical <- pData(gse_data) %>%
  dplyr::select(geo_accession, title) %>%
  mutate(
    # CRITICAL FIX: Decoded the author's cryptic _p (pre) and _d (during) labels
    Treatment_State = case_when(
      str_detect(title, "_p") ~ "Pre-Treatment",
      str_detect(title, "_d") ~ "Mid-Irradiation",
      TRUE ~ "Unknown"
    ),
    Patient_ID = str_extract(title, "^[0-9]+") 
  )

cat("--- 3. Processing & Translating Expression Data ---\n")
expr_matrix <- exprs(gse_data)
feature_data <- fData(gse_data)

if(max(expr_matrix, na.rm = TRUE) > 50) {
  cat("Applying log2(x+1) transform...\n")
  expr_matrix <- log2(expr_matrix + 1)
}

cat("Translating GenBank IDs (GB_LIST) to modern Gene Symbols...\n")
gb_ids <- as.character(feature_data$GB_LIST)
gb_ids_clean <- sapply(strsplit(gb_ids, ","), `[`, 1)

gb_ids_clean <- sub("\\..*", "", gb_ids_clean)
gb_ids_clean <- trimws(gb_ids_clean)

symbols <- rep(NA_character_, length(gb_ids_clean))
valid_gb <- !is.na(gb_ids_clean) & gb_ids_clean != ""

symbols[valid_gb] <- mapIds(
  org.Hs.eg.db,
  keys = gb_ids_clean[valid_gb],
  column = "SYMBOL",
  keytype = "ACCNUM",
  multiVals = "first"
)

valid_idx <- which(!is.na(symbols))
expr_mapped <- expr_matrix[valid_idx, , drop = FALSE]
mapped_symbols <- symbols[valid_idx]

cat("Aggregating translated probes (with NA protection)...\n")
expr_df <- as.data.frame(expr_mapped)
expr_df$Symbol <- mapped_symbols

expr_agg <- expr_df %>%
  group_by(Symbol) %>%
  summarise(across(everything(), ~ suppressWarnings(max(.x, na.rm = TRUE)))) %>%
  ungroup()

expr_agg_mat <- as.matrix(expr_agg[, -1])
expr_agg_mat[is.infinite(expr_agg_mat)] <- 0

expr_t <- as.data.frame(t(expr_agg_mat))
colnames(expr_t) <- expr_agg$Symbol
expr_t$geo_accession <- rownames(expr_t)

cohort <- merge(clinical, expr_t, by = "geo_accession")
cohort$OncotreeLineage <- "Cervix"
cohort$ModelID <- as.character(cohort$geo_accession)

cat("--- 4. Direct Type-Safe Data Injection ---\n")
current_cols <- colnames(cohort)
clean_current_cols <- gsub("[^A-Z0-9]", "", toupper(current_cols))
matched_count <- 0

for (expected_col in expected_vars) {
  if(expected_col %in% c("ModelID", "OncotreeLineage", "PrimaryOrMetastasis", "PRKDC_effect")) next
  if (endsWith(tolower(expected_col), "_cn") || endsWith(tolower(expected_col), "_mut")) next
  
  base_expected <- sub(" \\([0-9]+\\).*|\\.\\.[0-9]+\\..*", "", expected_col)
  clean_expected <- gsub("[^A-Z0-9]", "", toupper(base_expected))
  match_idx <- which(clean_current_cols == clean_expected)
  
  if (length(match_idx) > 0) {
    colnames(cohort)[match_idx[1]] <- expected_col
    clean_current_cols[match_idx[1]] <- "MAPPED" 
    matched_count <- matched_count + 1
  }
}
cat(sprintf("Successfully translated and mapped %d biological features!\n", matched_count))

cat("Constructing perfectly typed dataframe...\n")
n_patients <- nrow(cohort)
aligned_list <- list()

for (col in expected_vars) {
  if (col %in% colnames(cohort)) {
    val <- cohort[[col]]
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

final_df <- as_tibble(aligned_list)
final_df$geo_accession <- cohort$geo_accession
final_df$Treatment_State <- cohort$Treatment_State
final_df$Patient_ID <- cohort$Patient_ID

cat("--- 5. Prediction & Paired Dynamics Analysis ---\n")
var_check <- apply(final_df %>% dplyr::select(where(is.numeric)), 2, var, na.rm = TRUE)
active_features <- sum(var_check > 0, na.rm = TRUE)
cat(sprintf("Matrix check: %d features have active biological variance.\n", active_features))

preds <- predict(locked_model, new_data = final_df)

results <- final_df %>%
  dplyr::select(geo_accession, Patient_ID, Treatment_State) %>%
  mutate(Predicted_Dependency_Score = preds[[1]])

pred_sd <- sd(results$Predicted_Dependency_Score, na.rm = TRUE)
cat(sprintf("Prediction Variance (SD): %.5f\n", pred_sd))

write_csv(results, "data/processed/05_gse3578_predictions.csv")
cat("\n=============================================================\n")
cat("SUCCESS: Saved GSE3578 Paired Dynamics predictions to:\n")
cat("'data/processed/05_gse3578_predictions.csv'\n")
cat("=============================================================\n")