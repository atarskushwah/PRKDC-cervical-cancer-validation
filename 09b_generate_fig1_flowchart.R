# ==============================================================================
# scripts/09b_generate_fig1_flowchart.R
# ==============================================================================
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

cat("Generating clean editorial flowchart for Figure 1...\n")

flowchart <- grViz("
digraph study_flow {
  # Increased ranksep and nodesep for better spacing
  graph [layout = dot, rankdir = TB, nodesep = 0.5, ranksep = 0.7, fontname = 'Helvetica-Bold', bgcolor = 'white']
  
  node [shape = rectangle, style = 'filled,rounded', fontname = 'Helvetica', fontsize = 10, margin = '0.2,0.15', penwidth = 1.2]
  edge [color = '#333333', penwidth = 1.2, arrowhead = vee, arrowsize = 0.8]

  subgraph cluster_dev {
    label = 'I. MODEL DEVELOPMENT & LOCKING (DepMap)'
    labelloc = 't'
    style = 'dotted'; color = '#2B4C7E'; fontname = 'Helvetica-Bold'; fontsize = 11

    A [label = 'DepMap Discovery Cohort\n• Genome-scale CRISPR PRKDC Dependency\n• Matched Transcriptomics (N = 1,207)', fillcolor = '#F0F4F8', color = '#2B4C7E']
    B [label = 'Elastic-Net Penalized Regression\n• L1/L2 Regularization (Penalty = 0.01)\n• Cross-validated R² = 0.0190 | RMSE = 0.1301', fillcolor = '#FFF9E6', color = '#B78103']
    C [label = 'Locked PRKDC Dependency Signature\n• 19,074 Mapped Transcriptomic Features', fillcolor = '#EBF6ED', color = '#1E7E34']
    
    A -> B -> C
  }

  subgraph cluster_valid {
    label = 'II. CLINICAL TRANSLATION & VALIDATION'
    labelloc = 'b'  # <--- FIX: This moves the label to the bottom, avoiding the arrows
    style = 'dotted'; color = '#4A154B'; fontname = 'Helvetica-Bold'; fontsize = 11

    D [label = 'TCGA-CESC Primary\n• N = 309 Tumors\n• Purity Independent\n  (R = -0.2271, p = 6.11e-05)', fillcolor = '#F9F0F8', color = '#4A154B']
    E [label = 'GSE56363 Response\n• N = 21 Patients\n• CRT Outcome Metric\n  (OR = 1.64, AUC = 0.62)', fillcolor = '#F9F0F8', color = '#4A154B']
    F [label = 'GSE3578 Dynamics\n• N = 39 Paired Biopsies\n• Pre vs. Mid-Irradiation\n  (Wilcoxon p = 0.387)', fillcolor = '#F9F0F8', color = '#4A154B']
    G [label = 'GSE224327 Single-Cell\n• N = 12,631 Malignant Cells\n• Stroma/Immune Excluded\n  (SD = 0.00872)', fillcolor = '#F9F0F8', color = '#4A154B']
  }

  C -> D; C -> E; C -> F; C -> G
}
")

# Export to SVG -> PNG -> Convert to high-res TIFF natively on Mac -> Delete temporary PNG
svg_code <- export_svg(flowchart)
tmp_png <- "manuscript_figures/temp_fig1.png"
final_tiff <- "manuscript_figures/Fig1_Study_Flow.tiff"

rsvg_png(charToRaw(svg_code), file = tmp_png, width = 3000, height = 2400)
system(sprintf("sips -s format tiff -s dpiHeight 300 -s dpiWidth 300 %s --out %s", tmp_png, final_tiff), ignore.stdout = TRUE)
unlink(tmp_png)

cat("SUCCESS: Fixed high-res TIFF flowchart saved to 'manuscript_figures/Fig1_Study_Flow.tiff'!\n")