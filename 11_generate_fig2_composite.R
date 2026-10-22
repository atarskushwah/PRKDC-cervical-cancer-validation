# ==============================================================================
# scripts/11_generate_fig2_composite.R
# ==============================================================================
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)

setwd("~/CERV-DPK_Multiomics")

cat("Loading 560MB DepMap Discovery dataset...\n")
df <- read.csv("data/processed/01_harmonized_depmap_discovery.csv")

# Clean out any empty rows to prevent NA factor errors
df <- df %>% filter(!is.na(OncotreeLineage), OncotreeLineage != "")

# ------------------------------------------------------------------------------
# PANEL A: Pan-cancer landscape by lineage
# ------------------------------------------------------------------------------
lineage_medians <- df %>%
  group_by(OncotreeLineage) %>%
  summarize(median_dep = median(PRKDC_effect, na.rm = TRUE)) %>%
  arrange(median_dep)

df$OncotreeLineage <- factor(df$OncotreeLineage, levels = unique(lineage_medians$OncotreeLineage))

pA <- ggplot(df, aes(x = PRKDC_effect, y = OncotreeLineage)) +
  geom_boxplot(aes(fill = ifelse(OncotreeLineage == "Cervix", "Cervix", "Other")), 
               outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Cervix" = "#D95F02", "Other" = "#E0E0E0")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.y = element_text(size = 8)) +
  labs(x = "PRKDC Dependency (CRISPR Gene Effect)", 
       y = "Cancer Lineage", title = "Pan-cancer PRKDC Dependency")

# ------------------------------------------------------------------------------
# PANEL B: Cervical cell-line positions
# ------------------------------------------------------------------------------
cervix_df <- df %>%
  filter(OncotreeLineage == "Cervix") %>%
  group_by(CellLineName) %>%
  summarize(PRKDC_effect = mean(PRKDC_effect, na.rm = TRUE)) %>%
  arrange(PRKDC_effect)

# The levels are now guaranteed to be unique
cervix_df$CellLineName <- factor(cervix_df$CellLineName, levels = unique(cervix_df$CellLineName))

pB <- ggplot(cervix_df, aes(x = PRKDC_effect, y = CellLineName)) +
  geom_segment(aes(x = 0, xend = PRKDC_effect, y = CellLineName, yend = CellLineName), color = "grey50") +
  geom_point(color = "#D95F02", size = 3) +
  theme_classic(base_size = 12) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "CRISPR Gene Effect", y = "", title = "Cervical Cell Lines")

# ------------------------------------------------------------------------------
# PANEL C: Expression vs. Gene Effect
# ------------------------------------------------------------------------------
pC <- ggplot(df, aes(x = PRKDC_expression, y = PRKDC_effect)) +
  geom_point(alpha = 0.3, color = "#2B4C7E") +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  theme_classic(base_size = 12) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  labs(x = "PRKDC mRNA Expression", y = "CRISPR Gene Effect", title = "Expression vs. Dependency")

# ------------------------------------------------------------------------------
# STITCH & EXPORT
# ------------------------------------------------------------------------------
cat("Stitching panels A, B, and C...\n")
composite_fig <- (pA) / (pB | pC) + 
  plot_annotation(tag_levels = 'A') + 
  plot_layout(heights = c(2, 1))

cat("Saving as high-resolution TIFF...\n")
dir.create("manuscript_figures", showWarnings = FALSE)
ggsave("manuscript_figures/Fig2_PanCancer_Composite.tiff", 
       plot = composite_fig, device = "tiff", 
       width = 10, height = 12, dpi = 300, compression = "lzw")

cat("SUCCESS: Figure 2 saved to 'manuscript_figures/Fig2_PanCancer_Composite.tiff'!\n")