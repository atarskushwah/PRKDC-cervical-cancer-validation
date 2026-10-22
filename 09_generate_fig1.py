# ==============================================================================
# scripts/09_generate_fig1.py
# ==============================================================================
import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.shapes import MSO_SHAPE
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

# Create output folder if it doesn't exist
os.makedirs("manuscript_figures", exist_ok=True)

prs = Presentation()
prs.slide_width = Inches(11)
prs.slide_height = Inches(8.5)
slide = prs.slides.add_slide(prs.slide_layouts[6])  # Blank layout

def add_box(slide, left, top, width, height, text, bg_rgb, border_rgb):
    shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = bg_rgb
    shape.line.color.rgb = border_rgb
    shape.line.width = Pt(1.5)
    
    tf = shape.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(11)
    p.font.color.rgb = RGBColor(0, 0, 0)
    p.alignment = PP_ALIGN.CENTER
    return shape

# --- Top Section: Phases 1 & 2 (Model Development) ---
add_box(slide, Inches(3.5), Inches(0.5), Inches(4.0), Inches(0.9), 
        "DepMap Discovery Cohort\nCRISPR PRKDC Dependency (N = 1,207)", 
        RGBColor(227, 242, 253), RGBColor(21, 101, 192))

add_box(slide, Inches(3.5), Inches(1.8), Inches(4.0), Inches(0.9), 
        "Elastic-Net Penalized Regression\n(Penalty = 0.01 | CV-R² = 0.0190)", 
        RGBColor(255, 249, 196), RGBColor(251, 192, 45))

add_box(slide, Inches(3.5), Inches(3.1), Inches(4.0), Inches(0.9), 
        "Locked PRKDC Dependency Signature\n(19,074 Mapped Transcripts)", 
        RGBColor(200, 230, 201), RGBColor(46, 125, 50))

# --- Bottom Section: Phases 3 & 4 (Clinical Validation) ---
cohorts = [
    ("TCGA-CESC Primary", "N = 309 Tumors\nCPE Check (R = -0.2271)"),
    ("GSE56363 CRT Response", "N = 21 Patients\nOR = 1.64, AUC = 0.62"),
    ("GSE3578 Dynamics", "N = 39 Paired Biopsies\nPre vs. Mid-Irradiation"),
    ("GSE224327 Single-Cell", "N = 12,631 Malignant Cells\n(Stroma/Immune Excluded)")
]

start_left = Inches(0.4)
spacing = Inches(2.6)

for i, (title, desc) in enumerate(cohorts):
    add_box(
        slide, 
        start_left + (i * spacing), 
        Inches(4.8), 
        Inches(2.4), 
        Inches(1.5), 
        f"{title}\n\n{desc}", 
        RGBColor(243, 229, 245), 
        RGBColor(123, 31, 162)
    )

prs.save("manuscript_figures/Fig1_Study_Flow.pptx")
print("SUCCESS: Editable PowerPoint saved to 'manuscript_figures/Fig1_Study_Flow.pptx'")
