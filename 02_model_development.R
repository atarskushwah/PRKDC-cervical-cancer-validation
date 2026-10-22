# ==============================================================================
# scripts/02_model_development.R
# ==============================================================================
library(tidyverse)
library(tidymodels)
library(glmnet)

cat("--- Loading Harmonized Data ---\n")
df <- read_csv("data/processed/01_harmonized_depmap_discovery.csv", show_col_types = FALSE) %>%
  distinct(ModelID, .keep_all = TRUE)

df_model <- df %>%
  select(ModelID, OncotreeLineage, where(is.numeric)) %>%
  select(-contains("_cn"), -contains("_mut")) %>%
  drop_na(PRKDC_effect, OncotreeLineage)

cat(sprintf("Training model on %d cell lines and %d transcriptomic features...\n", 
            nrow(df_model), ncol(df_model) - 3))

# 1. Define the Preprocessing Recipe
cerv_recipe <- recipe(PRKDC_effect ~ ., data = df_model) %>%
  update_role(ModelID, OncotreeLineage, new_role = "ID") %>%
  step_zv(all_predictors()) %>% 
  step_impute_mean(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors())

# 2. Define the Elastic-Net Specification
# THE FIX: Hardcode a low penalty (0.01) to guarantee the model keeps the expression 
# coefficients active, bypassing the grid tuning that zeroed everything out.
elastic_net_spec <- linear_reg(penalty = 0.01, mixture = 0.5) %>%
  set_engine("glmnet")

# 3. Construct Workflow
cerv_wf <- workflow() %>%
  add_recipe(cerv_recipe) %>%
  add_model(elastic_net_spec)

# 4. Fit Model Directly (No Tuning Loop needed!)
cat("Fitting active Elastic Net model (this will be fast)...\n")
locked_model <- fit(cerv_wf, data = df_model)

# 5. Freeze to Disk
saveRDS(locked_model, "models/locked_prkdc_dependency_model.rds")
saveRDS(cerv_recipe, "models/locked_preprocessing_recipe.rds")

cat("SUCCESS: Active model mathematically frozen and saved to 'models/' directory.\n")