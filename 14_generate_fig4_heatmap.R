library(ComplexHeatmap)
library(circlize)
library(dplyr)

setwd("~/CERV-DPK_Multiomics")

expr_file <- "data/processed/tcga_cesc_expression.rds"
clin_file <- "data/processed/03_tcga_cesc_predictions.csv"

if(!file.exists(expr_file)) stop("Missing TCGA expression data.")
if(!file.exists(clin_file)) stop("Missing clinical predictions data.")

tcga_expr <- readRDS(expr_file)
tcga_clin <- read.csv(clin_file)

colnames(tcga_expr) <- substr(colnames(tcga_expr), 1, 12)
rownames(tcga_clin) <- substr(tcga_clin$patient, 1, 12)
common <- intersect(colnames(tcga_expr), rownames(tcga_clin))

tcga_expr <- tcga_expr[, common]
tcga_clin <- tcga_clin[common, ]
tcga_clin <- tcga_clin[order(tcga_clin$Predicted_Dependency_Score), ]
tcga_expr <- tcga_expr[, rownames(tcga_clin)]

top_genes <- names(sort(apply(tcga_expr, 1, var, na.rm = TRUE), decreasing = TRUE)[1:1500])
mat_scaled <- t(scale(t(tcga_expr[top_genes, ])))

n_pat <- nrow(tcga_clin)

get_safe_breaks <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if(rng[1] == rng[2]) rng <- c(rng[1] - 1, rng[1] + 1)
  return(rng)
}

dep_bks <- get_safe_breaks(tcga_clin$Predicted_Dependency_Score)

# Pull real attributes from your metadata frame rather than generating samples
ha <- HeatmapAnnotation(
  Dependency = tcga_clin$Predicted_Dependency_Score,
  HR_Pathway = tcga_clin$HR_Pathway,
  NHEJ_Pathway = tcga_clin$NHEJ_Pathway,
  Checkpoint = tcga_clin$Checkpoint,
  Rep_Stress = tcga_clin$Rep_Stress,
  HPV_Group = tcga_clin$HPV_Group,
  Histology = tcga_clin$Histology,
  Methylation = tcga_clin$Methylation_Status,
  col = list(
    Dependency = colorRamp2(dep_bks, c("white", "#D95F02")),
    HR_Pathway = colorRamp2(c(-2, 0, 2), c("#2B4C7E", "white", "#C62828")),
    NHEJ_Pathway = colorRamp2(c(-2, 0, 2), c("#2B4C7E", "white", "#C62828")),
    Checkpoint = colorRamp2(c(-2, 0, 2), c("#2B4C7E", "white", "#C62828")),
    Rep_Stress = colorRamp2(c(-2, 0, 2), c("#2B4C7E", "white", "#C62828")),
    HPV_Group = c("HPV16" = "#1b9e77", "HPV18" = "#d95f02", "Other HPV" = "#7570b3", "Negative" = "#e7298a"),
    Histology = c("Squamous" = "#e6ab02", "Adenocarcinoma" = "#a6761d"),
    Methylation = c("CIMP-High" = "#66a61e", "CIMP-Low" = "#e41a1c")
  ),
  annotation_name_side = "left"
)

tiff("manuscript_figures/Fig4_RepairState_Atlas.tiff", width = 4200, height = 3600, res = 300, compression = "lzw")
ht <- Heatmap(mat_scaled, 
              name = "Z-score",
              top_annotation = ha,
              show_row_names = FALSE, 
              show_column_names = FALSE,
              cluster_columns = FALSE, 
              cluster_rows = TRUE,
              use_raster = TRUE,
              col = colorRamp2(c(-2, 0, 2), c("#1565C0", "white", "#C62828")))
draw(ht)
dev.off()
cat("SUCCESS: Comprehensive Figure 4 saved safely!\n")