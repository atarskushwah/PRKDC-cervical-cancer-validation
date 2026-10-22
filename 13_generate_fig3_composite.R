# ==============================================================================
# scripts/13_generate_fig3_composite.R
# ==============================================================================
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(tidymodels)
library(broom)

setwd("~/CERV-DPK_Multiomics")

cat("Loading discovery dataset and locked model...\n")
df <- read.csv("data/processed/01_harmonized_depmap_discovery.csv")
model <- readRDS("models/locked_prkdc_dependency_model.rds")

# ------------------------------------------------------------------------------
# PANEL 3A: Nested Cross-Validation Performance
# ------------------------------------------------------------------------------
# Extract true CV metrics from the fitted workflow resamples
cv_metrics <- collect_metrics(model)
cv_data <- cv_metrics %>% filter(.metric == "rsq") %>% rename(Rsquared = mean)
cv_data$Resample <- paste0("Fold", seq_len(nrow(cv_data)))

pA <- ggplot(cv_data, aes(y = Rsquared, x = "Outer CV Folds")) +
  geom_boxplot(fill = "#EBF6ED", color = "#1E7E34", width = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2, color = "black", alpha = 0.6) +
  theme_classic(base_size = 12) +
  labs(x = "", y = bquote("Cross-validated "~R^2), title = "Model Robustness")

# ------------------------------------------------------------------------------
# PANEL 3B: Observed vs. Predicted Dependency
# ------------------------------------------------------------------------------
# Extract predictions strictly using new_data
pred_vals <- predict(model, new_data = df)$.pred
df$Predicted_Score <- pred_vals

pB <- ggplot(df, aes(x = PRKDC_effect, y = Predicted_Score)) +
  geom_point(alpha = 0.4, color = "#2B4C7E") +
  geom_smooth(method = "lm", color = "black", linetype = "solid", linewidth = 1) +
  theme_classic(base_size = 12) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
  labs(x = "Observed PRKDC Effect (CRISPR)", y = "Predicted Dependency Score", title = "Prediction Accuracy")

# ------------------------------------------------------------------------------
# PANEL 3C: Stable Feature Coefficients (Top 15)
# ------------------------------------------------------------------------------
# Extract active elastic-net coefficients rather than using varImp
coef_df <- broom::tidy(extract_fit_parsnip(model)) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  arrange(desc(abs(estimate))) %>%
  head(15)

pC <- ggplot(coef_df, aes(x = reorder(term, abs(estimate)), y = abs(estimate))) +
  geom_bar(stat = "identity", fill = "#B78103", alpha = 0.8) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = "Transcriptomic Feature", y = "Absolute Coefficient", title = "Feature Stability")

# ------------------------------------------------------------------------------
# STITCH & EXPORT
# ------------------------------------------------------------------------------
cat("Stitching composite figure...\n")
composite_fig3 <- (pA | pB) / pC + plot_annotation(tag_levels = 'A')

dir.create("manuscript_figures", showWarnings = FALSE)
ggsave("manuscript_figures/Fig3_Model_Composite.tiff", 
       plot = composite_fig3, device = "tiff", 
       width = 12, height = 8, dpi = 300, compression = "lzw")

cat("SUCCESS: Figure 3 saved to 'manuscript_figures/Fig3_Model_Composite.tiff'!\n")