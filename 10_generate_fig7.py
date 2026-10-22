# ==============================================================================
# scripts/10_generate_fig7.py
# ==============================================================================
import os
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.enum.dml import MSO_LINE

os.makedirs("manuscript_figures", exist_ok=True)
prs = Presentation()
prs.slide_width = Inches(10)
prs.slide_height = Inches(7.5)
slide = prs.slides.add_slide(prs.slide_layouts[6])

def add_node(slide, left, top, width, height, text, bg_rgb, border_rgb):
    shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = bg_rgb
    shape.line.color.rgb = border_rgb
    shape.line.width = Pt(1.5)
    tf = shape.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(12)
    p.font.color.rgb = RGBColor(0, 0, 0)
    p.alignment = PP_ALIGN.CENTER
    return shape

def add_arrow(slide, start_x, start_y, end_x, end_y, dashed=False):
    # FIX: Use MSO_CONNECTOR.STRAIGHT for connecting lines
    connector = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, start_x, start_y, end_x, end_y)
    connector.line.color.rgb = RGBColor(50, 50, 50)
    connector.line.width = Pt(2.5)
    if dashed:
        connector.line.dash_style = MSO_LINE.DASH
    return connector

# --- Path 1: The Hypothesis (HR Deficiency) ---
add_node(slide, Inches(1.5), Inches(1.0), Inches(2.5), Inches(0.8), "Oncogenic Replication Stress\n/ Chromatin Alterations", RGBColor(255, 235, 238), RGBColor(198, 40, 40))
add_node(slide, Inches(1.5), Inches(3.0), Inches(2.5), Inches(0.8), "Homologous Recombination\n(HR) Impairment", RGBColor(255, 235, 238), RGBColor(198, 40, 40))

# Dashed lines for hypothesis
add_arrow(slide, Inches(2.75), Inches(1.8), Inches(2.75), Inches(3.0), dashed=True)
add_arrow(slide, Inches(2.75), Inches(3.8), Inches(4.5), Inches(5.0), dashed=True)

# --- Path 2: The Data-Supported Mechanism (NHEJ Dependency) ---
add_node(slide, Inches(6.0), Inches(1.0), Inches(2.5), Inches(0.8), "DNA Double-Strand Breaks\n(Exacerbated by Radiation)", RGBColor(227, 242, 253), RGBColor(21, 101, 192))
add_node(slide, Inches(6.0), Inches(3.0), Inches(2.5), Inches(0.8), "PRKDC (DNA-PK) Upregulation\n& Assembly at Break Sites", RGBColor(232, 245, 233), RGBColor(46, 125, 50))
add_node(slide, Inches(3.75), Inches(5.0), Inches(2.5), Inches(0.8), "Absolute Reliance on\nDNA-PK-mediated NHEJ", RGBColor(255, 243, 224), RGBColor(239, 108, 0))

# Solid lines for data-backed findings
add_arrow(slide, Inches(7.25), Inches(1.8), Inches(7.25), Inches(3.0), dashed=False)
add_arrow(slide, Inches(7.25), Inches(3.8), Inches(5.5), Inches(5.0), dashed=False)

# Legend Box
leg = add_node(slide, Inches(7.0), Inches(6.0), Inches(2.5), Inches(1.0), "Legend:\n▬ Data-supported link\n--- Proposed mechanism", RGBColor(245, 245, 245), RGBColor(150, 150, 150))
leg.text_frame.paragraphs[0].font.size = Pt(10)

prs.save("manuscript_figures/Fig7_Biological_Model.pptx")
print("SUCCESS: Editable Figure 7 saved to 'manuscript_figures/Fig7_Biological_Model.pptx'")
