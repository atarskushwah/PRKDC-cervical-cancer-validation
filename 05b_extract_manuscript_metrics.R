# ==============================================================================
# scripts/05b_extract_manuscript_metrics.R
# ==============================================================================
library(tidyverse)
library(tidymodels)

cat("=============================================================\n")
cat("--- 1. Extracting Training Metrics (R-squared & RMSE) ---\n")
cat("=============================================================\n")

# Load harmonized training data
df <- read_csv("data/processed/01_harmonized_depmap_discovery.csv", show_col_types = FALSE) %>%
  distinct(ModelID, .keep_all = TRUE) %>%
  dplyr::select(ModelID, OncotreeLineage, where(is.numeric)) %>% # FIXED: Explicit dplyr namespace
  dplyr::select(-contains("_cn"), -contains("_mut")) %>%         # FIXED: Explicit dplyr namespace
  drop_na(PRKDC_effect, OncotreeLineage)

# Load the locked recipe
locked_recipe <- readRDS("models/locked_preprocessing_recipe.rds")

# Re-declare the locked model parameters (Penalty = 0.01)
locked_spec <- linear_reg(penalty = 0.01, mixture = 0.5) %>%
  set_engine("glmnet")

# 5-fold Cross Validation
set.seed(2026)
cerv_folds <- vfold_cv(df, v = 5, strata = OncotreeLineage)

wf <- workflow() %>%
  add_recipe(locked_recipe) %>%
  add_model(locked_spec)

cat("Running 5-fold CV to calculate formal error metrics (this will take ~30 seconds)...\n")
cv_res <- fit_resamples(wf, resamples = cerv_folds)
metrics <- collect_metrics(cv_res)

rmse_val <- metrics %>% filter(.metric == "rmse") %>% pull(mean)
r2_val <- metrics %>% filter(.metric == "rsq") %>% pull(mean)

cat(sprintf("-> Cross-Validated RMSE: %.4f\n", rmse_val))
cat(sprintf("-> Cross-Validated R-squared: %.4f\n", r2_val))


cat("\n=============================================================\n")
cat("--- 2. Extracting GSE3578 Paired Shift Statistics ---\n")
cat("=============================================================\n")

# Load the clinical predictions we just generated
gse3578 <- read_csv("data/processed/05_gse3578_predictions.csv", show_col_types = FALSE)

# Filter to valid states and ensure we only keep patients with strictly paired biopsies
paired_data <- gse3578 %>%
  filter(Treatment_State %in% c("Pre-Treatment", "Mid-Irradiation")) %>%
  group_by(Patient_ID) %>%
  filter(n_distinct(Treatment_State) == 2) %>%
  ungroup()

# Pivot data wide to align the pairs row-by-row
paired_wide <- paired_data %>%
  dplyr::select(Patient_ID, Treatment_State, Predicted_Dependency_Score) %>% # FIXED: Explicit dplyr namespace
  pivot_wider(names_from = Treatment_State, values_from = Predicted_Dependency_Score, values_fn = mean)

cat(sprintf("Found %d perfectly paired patient biopsies.\n", nrow(paired_wide)))

# Calculate means
mean_pre <- mean(paired_wide$`Pre-Treatment`, na.rm = TRUE)
mean_mid <- mean(paired_wide$`Mid-Irradiation`, na.rm = TRUE)

cat(sprintf("-> Mean Pre-Treatment Score:  %.5f\n", mean_pre))
cat(sprintf("-> Mean Mid-Irradiation Score: %.5f\n", mean_mid))

# Run paired Wilcoxon signed-rank test
res <- wilcox.test(paired_wide$`Mid-Irradiation`, paired_wide$`Pre-Treatment`, paired = TRUE)

cat(sprintf("-> Paired Wilcoxon p-value: %.2e\n", res$p.value))
cat("=============================================================\n")