# ==============================================================================
# scripts/08_generate_figures.R
# ==============================================================================
# Run this in your console if you are missing any packages: 
# install.packages(c("patchwork", "pROC", "ggpubr"))

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(pROC)

dir.create("manuscript_figures", showWarnings = FALSE)

cat("Generating Figure 2: Pan-Cancer Landscape (DepMap)...\n")
# ------------------------------------------------------------------------------
depmap <- read_csv("data/processed/01_harmonized_depmap_discovery.csv", show_col_types = FALSE) %>% 
  drop_na(PRKDC_effect)

# Calculate median dependency per lineage to sort the plot nicely
depmap_summary <- depmap %>%
  group_by(OncotreeLineage) %>%
  mutate(Median_Effect = median(PRKDC_effect, na.rm = TRUE)) %>%
  ungroup()

fig2 <- ggplot(depmap_summary, aes(x = reorder(OncotreeLineage, Median_Effect), y = PRKDC_effect)) +
  geom_boxplot(aes(fill = OncotreeLineage == "Cervix"), outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("gray80", "#D55E00")) +
  coord_flip() +
  labs(x = "Cancer Lineage", y = "PRKDC Dependency (CRISPR Gene Effect)", 
       title = "Pan-cancer PRKDC Dependency") +
  theme_classic() +
  theme(legend.position = "none", axis.text.y = element_text(size = 6),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("manuscript_figures/Fig2_PanCancer_Landscape.tiff", plot = fig2, width = 7, height = 8, dpi = 300, compression = "lzw")


cat("Generating Figure 3: Model Performance (DepMap)...\n")
# ------------------------------------------------------------------------------
locked_model <- readRDS("models/locked_prkdc_dependency_model.rds")
depmap_preds <- predict(locked_model, new_data = depmap)
depmap_plot <- depmap %>% mutate(Predicted = depmap_preds$.pred)

fig3 <- ggscatter(depmap_plot, x = "PRKDC_effect", y = "Predicted",
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
                  xlab = "Observed PRKDC Effect (CRISPR)", ylab = "Predicted Score",
                  title = "Observed vs. Predicted Dependency", color = "#0073C2FF", alpha = 0.4) + theme_classic()

ggsave("manuscript_figures/Fig3_Model_Performance.tiff", plot = fig3, width = 6, height = 5, dpi = 300, compression = "lzw")


cat("Generating Figure 5: Clinical Response & Dynamics...\n")
# ------------------------------------------------------------------------------
# 5A: GSE56363 Response
gse56363 <- read_csv("data/processed/04_gse56363_predictions.csv", show_col_types = FALSE)
resp_col <- grep("response|Outcome", colnames(gse56363), ignore.case = TRUE, value = TRUE)[1]
gse_clean <- gse56363 %>% drop_na(!!sym(resp_col))
fig5a <- ggboxplot(gse_clean, x = resp_col, y = "Predicted_Dependency_Score", fill = resp_col, 
                   palette = c("#E69F00", "#56B4E9"), xlab = "Clinical Response", ylab = "Locked Score") +
  stat_compare_means() + theme_classic() + theme(legend.position = "none")

# 5B: GSE3578 Paired Shifts
gse3578 <- read_csv("data/processed/05_gse3578_predictions.csv", show_col_types = FALSE) %>%
  filter(Treatment_State %in% c("Pre-Treatment", "Mid-Irradiation"))
gse3578$Treatment_State <- factor(gse3578$Treatment_State, levels = c("Pre-Treatment", "Mid-Irradiation"))
fig5b <- ggpaired(gse3578, x = "Treatment_State", y = "Predicted_Dependency_Score", id = "Patient_ID", 
                  color = "Treatment_State", line.color = "gray", palette = c("#0073C2FF", "#EFC000FF"),
                  xlab = "Biopsy State", ylab = "Predicted Score") +
  stat_compare_means(paired = TRUE) + theme_classic() + theme(legend.position = "none")

fig5_combined <- fig5a + fig5b + plot_annotation(tag_levels = 'A')
ggsave("manuscript_figures/Fig5_Clinical_Dynamics.tiff", plot = fig5_combined, width = 10, height = 5, dpi = 300, compression = "lzw")


cat("Generating Figure 6: Single-Cell Malignant Localization...\n")
# ------------------------------------------------------------------------------
gse224327 <- read_csv("data/processed/06_gse224327_predictions.csv", show_col_types = FALSE)
fig6 <- ggplot(gse224327, aes(x = Patient_ID, y = Predicted_Dependency_Score, fill = Patient_ID)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_viridis_d(option = "plasma") + labs(title = "Malignant-Cell Distribution", x = "Patient", y = "Score") +
  theme_classic() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("manuscript_figures/Fig6_SingleCell_Violins.tiff", plot = fig6, width = 7, height = 5, dpi = 300, compression = "lzw")

cat("SUCCESS: All 4 data-driven figures saved as high-res TIFFs to the 'manuscript_figures' folder!\n")