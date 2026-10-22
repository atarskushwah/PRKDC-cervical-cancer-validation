# ==============================================================================
# scripts/06b_gse224327_single_cell_localization.R
# ==============================================================================
library(Seurat)
library(tidyverse)
library(tidymodels)
library(Matrix) # CRITICAL: Used for memory-safe sparse transposing

cat("--- 1. Data & Model Initialization ---\n")
locked_recipe <- readRDS("models/locked_preprocessing_recipe.rds")
locked_model <- readRDS("models/locked_prkdc_dependency_model.rds")

expected_vars <- locked_recipe$var_info$variable
template_data <- locked_recipe$template

sc_dir <- "data/raw/GSE224327_scRNA"
files <- list.files(sc_dir, full.names = TRUE)
prefixes <- unique(sub("_barcodes.*|_features.*|_matrix.*", "", files))

epithelial_markers <- c("EPCAM", "KRT18", "KRT19", "CDH1")
stromal_immune_markers <- c("PTPRC", "CD3E", "COL1A1", "DCN")

results_list <- list()

cat("--- 2. Ultra-Low-Memory Patient Processing ---\n")

for(pref in prefixes) {
  pt_name <- str_extract(basename(pref), "PT[0-9]+")
  cat(sprintf("\n>>> Booting up %s...\n", pt_name))
  
  mat <- ReadMtx(
    mtx = paste0(pref, "_matrix.mtx.gz"),
    features = paste0(pref, "_features.tsv.gz"),
    cells = paste0(pref, "_barcodes.tsv.gz")
  )
  
  sc <- CreateSeuratObject(counts = mat, project = pt_name, min.cells = 3, min.features = 200)
  rm(mat) 
  invisible(gc())
  
  sc <- NormalizeData(sc, verbose = FALSE)
  sc <- AddModuleScore(sc, features = list(epithelial_markers), name = "Epithelial_Score")
  sc <- AddModuleScore(sc, features = list(stromal_immune_markers), name = "Stromal_Immune_Score")
  
  malignant <- subset(sc, Epithelial_Score1 > 0 & Stromal_Immune_Score1 < 0)
  n_malig <- ncol(malignant)
  cat(sprintf("  -> Isolated %d malignant cells.\n", n_malig))
  
  if(n_malig == 0) {
    rm(sc, malignant)
    invisible(gc())
    next
  }
  
  malignant <- JoinLayers(malignant)
  expr_matrix <- GetAssayData(malignant, layer = "data")
  
  shared_genes <- intersect(rownames(expr_matrix), sub(" \\([0-9]+\\).*|\\.\\.[0-9]+\\..*", "", expected_vars))
  
  # Step A: Filter while sparse
  sparse_subset <- expr_matrix[shared_genes, , drop = FALSE]
  
  # CRITICAL FIX: Transpose while sparse! Takes almost 0 RAM.
  sparse_subset_t <- Matrix::t(sparse_subset)
  
  # Instantly wipe everything massive from memory
  rm(sc, malignant, expr_matrix, sparse_subset)
  invisible(gc())
  
  # Step B: Inner Batching (1000 cells at a time)
  cell_barcodes <- rownames(sparse_subset_t)
  batch_size <- 1000
  n_batches <- ceiling(length(cell_barcodes) / batch_size)
  
  pt_results <- list()
  
  for(i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(cell_barcodes))
    bcs <- cell_barcodes[start_idx:end_idx]
    
    # Isolate chunk and coerce to dense dataframe (only ~10MB max!)
    chunk_sparse <- sparse_subset_t[bcs, , drop = FALSE]
    cohort <- as.data.frame(as.matrix(chunk_sparse))
    
    cohort$barcode <- rownames(cohort)
    cohort$Patient_ID <- pt_name
    cohort$OncotreeLineage <- "Cervix"
    cohort$ModelID <- cohort$barcode
    
    current_cols <- colnames(cohort)
    clean_current_cols <- gsub("[^A-Z0-9]", "", toupper(current_cols))
    
    for (expected_col in expected_vars) {
      if(expected_col %in% c("ModelID", "OncotreeLineage", "PrimaryOrMetastasis", "PRKDC_effect")) next
      if (endsWith(tolower(expected_col), "_cn") || endsWith(tolower(expected_col), "_mut")) next
      
      base_expected <- sub(" \\([0-9]+\\).*|\\.\\.[0-9]+\\..*", "", expected_col)
      clean_expected <- gsub("[^A-Z0-9]", "", toupper(base_expected))
      match_idx <- which(clean_current_cols == clean_expected)
      
      if (length(match_idx) > 0) {
        colnames(cohort)[match_idx[1]] <- expected_col
        clean_current_cols[match_idx[1]] <- "MAPPED" 
      }
    }
    
    n_cells_batch <- nrow(cohort)
    aligned_list <- list()
    
    for (col in expected_vars) {
      if (col %in% colnames(cohort)) {
        val <- cohort[[col]]
        if(is.data.frame(val) || is.matrix(val)) val <- val[[1]]
      } else {
        val <- rep(0, n_cells_batch)
      }
      
      if (!is.null(template_data) && col %in% colnames(template_data)) {
        tmpl_class <- class(template_data[[col]])[1]
        
        if (tmpl_class == "numeric") val <- as.numeric(val)
        else if (tmpl_class == "integer") val <- as.integer(val)
        else if (tmpl_class == "logical") val <- as.logical(val)
        else if (tmpl_class == "character") val <- as.character(val)
        else if (tmpl_class == "factor") val <- factor(as.character(val), levels = levels(template_data[[col]]))
        
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
    final_df$barcode <- cohort$barcode
    final_df$Patient_ID <- cohort$Patient_ID
    
    preds <- predict(locked_model, new_data = final_df)
    
    pt_results[[i]] <- final_df %>%
      dplyr::select(barcode, Patient_ID) %>%
      mutate(Predicted_Dependency_Score = preds[[1]])
    
    rm(chunk_sparse, cohort, aligned_list, final_df)
    invisible(gc())
  }
  
  results_list[[pt_name]] <- bind_rows(pt_results)
  cat(sprintf("  -> %s predictions complete. Wiping RAM...\n", pt_name))
  
  rm(sparse_subset_t, pt_results)
  invisible(gc())
}

cat("\n--- 3. Aggregation & Output ---\n")
final_results <- bind_rows(results_list)

pred_sd <- sd(final_results$Predicted_Dependency_Score, na.rm = TRUE)
cat(sprintf("Total Prediction Variance (SD): %.5f\n", pred_sd))

write_csv(final_results, "data/processed/06_gse224327_predictions.csv")
cat("\n=============================================================\n")
cat("SUCCESS: Saved Single-Cell Malignant predictions to:\n")
cat("'data/processed/06_gse224327_predictions.csv'\n")
cat("=============================================================\n")