library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("~/CERV-DPK_Multiomics")
sc_obj <- readRDS("data/processed/gse224327_seurat_object.rds")

if(any(duplicated(colnames(sc_obj)))) sc_obj <- RenameCells(sc_obj, make.unique = TRUE)

sc_obj <- NormalizeData(sc_obj, verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)

# Explicitly map patient identity from Seurat object to ensure discrete clusters
sc_obj$patient_id <- sc_obj$orig.ident

pA <- DimPlot(sc_obj, reduction = "umap", group.by = "compartment", pt.size = 0.4) +
  theme_classic(base_size = 11) +
  labs(title = "A. Major Compartments", x = "UMAP 1", y = "UMAP 2")

pB <- ggplot(sc_obj@meta.data, aes(x = compartment, y = Predicted_Dependency_Score, fill = compartment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none") +
  labs(title = "B. Score by Cell Type", x = "", y = "Locked Dependency Score")

pC <- DimPlot(sc_obj, reduction = "umap", group.by = "patient_id", pt.size = 0.3) +
  theme_classic(base_size = 11) +
  labs(title = "C. Patient Distribution", x = "UMAP 1", y = "UMAP 2")

# Calculate real pseudobulk metrics per patient
resp_df <- sc_obj@meta.data %>%
  group_by(patient_id, Clinical_Response) %>%
  summarize(Pseudobulk_Score = mean(Predicted_Dependency_Score, na.rm = TRUE), .groups = 'drop')

pD <- ggplot(resp_df, aes(x = Clinical_Response, y = Pseudobulk_Score, fill = Clinical_Response)) +
  geom_boxplot(width = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none") +
  labs(title = "D. Responder Pseudobulk", x = "Clinical Response", y = "Mean Pathway Score")

composite_fig6 <- (pA + pB) / (pC + pD)

tiff("manuscript_figures/Fig6_SingleCell_Composite.tiff", width = 4000, height = 3600, res = 300, compression = "lzw")
print(composite_fig6)
dev.off()
cat("SUCCESS: Complete 4-panel Figure 6 saved successfully!\n")