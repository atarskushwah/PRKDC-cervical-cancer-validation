# ==============================================================================
# scripts/12b_generate_fig5_complete.R
# ==============================================================================
library(ggplot2)
library(dplyr)
library(patchwork)
library(pROC)
library(ggpubr)
library(tidyr)

setwd("~/CERV-DPK_Multiomics")

cat("Loading complete clinical datasets...\n")
df_56363 <- read.csv("data/processed/04_gse56363_predictions.csv")
df_3578 <- read.csv("data/processed/05_gse3578_pathway_scores.csv")

# --- PANEL A: Pretreatment score ---
df_56363$response_label <- factor(df_56363$response_binary, 
                                  levels = c(0, 1), 
                                  labels = c("Non-Complete", "Complete"))

pA <- ggplot(df_56363, aes(x = response_label, y = Predicted_Dependency_Score, fill = response_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("Non-Complete" = "#56B4E9", "Complete" = "#E69F00")) +
  theme_classic(base_size = 12) + 
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test", label.x = 1.35, label.y = 0.62) +
  labs(x = "Clinical Response", y = "Locked Dependency Score", title = "GSE56363: CRT Response")

# --- PANEL B: ROC curve ---
roc_obj <- roc(df_56363$response_binary, df_56363$Predicted_Dependency_Score, quiet = TRUE)
ciobj <- ci.se(roc_obj, specificities = seq(0, 1, 0.05))
dat.ci <- data.frame(x = as.numeric(rownames(ciobj)), lower = ciobj[, 1], upper = ciobj[, 3])

pB <- ggroc(roc_obj, color = "#D95F02", size = 1.2, legacy.axes = TRUE) +
  geom_ribbon(data = dat.ci, aes(x = 1 - x, ymin = lower, ymax = upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  theme_classic(base_size = 12) +
  annotate("text", x = 0.75, y = 0.25, label = sprintf("AUC = %.2f", roc_obj$auc), size = 5) +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "Predictive Performance")

# --- PANEL C: Paired pre-/mid-irradiation ---
df_3578$Treatment_State <- factor(df_3578$Treatment_State, levels = c("Pre-Treatment", "Mid-Irradiation"))

pC <- ggplot(df_3578, aes(x = Treatment_State, y = Predicted_Dependency_Score, group = Patient_ID)) +
  geom_line(color = "grey70", alpha = 0.8) + 
  geom_point(aes(color = Treatment_State), size = 3) +
  scale_color_manual(values = c("Pre-Treatment" = "#1565C0", "Mid-Irradiation" = "#C62828")) +
  theme_classic(base_size = 12) + 
  theme(legend.position = "none") +
  annotate("text", x = 1.5, y = 0.24, label = "Wilcoxon, p = 0.387", size = 4.5) +
  labs(x = "Biopsy State", y = "Predicted Dependency Score", title = "GSE3578: Dependency Dynamics")

# --- PANEL D: Pathway-component dynamics ---
df_path <- df_3578 %>%
  select(Patient_ID, Treatment_State, HR_Score, NHEJ_Score) %>%
  pivot_longer(cols = c("HR_Score", "NHEJ_Score"), names_to = "Pathway", values_to = "Score")

df_path$Pathway <- factor(df_path$Pathway, levels = c("HR_Score", "NHEJ_Score"), 
                          labels = c("Homologous Recombination", "NHEJ (DNA-PK)"))

pD <- ggplot(df_path, aes(x = Pathway, y = Score, fill = Treatment_State)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  scale_fill_manual(values = c("Pre-Treatment" = "#1565C0", "Mid-Irradiation" = "#C62828")) +
  theme_classic(base_size = 12) + 
  theme(legend.position = "bottom", legend.title = element_blank()) +
  stat_compare_means(aes(group = Treatment_State), label = "p.signif", method = "wilcox.test") +
  labs(x = "DNA Repair Pathway", y = "Pathway Activity (Z-score)", title = "Pathway Component Dynamics")

# --- STITCH & EXPORT ---
cat("Stitching composite figure...\n")
composite_fig5 <- (pA | pB) / (pC | pD) + plot_annotation(tag_levels = 'A')

ggsave("manuscript_figures/Fig5_Clinical_Composite.tiff", 
       plot = composite_fig5, device = "tiff", 
       width = 12, height = 10, dpi = 300, compression = "lzw")

cat("SUCCESS: Corrected Figure 5 saved to 'manuscript_figures/Fig5_Clinical_Composite.tiff'!\n")