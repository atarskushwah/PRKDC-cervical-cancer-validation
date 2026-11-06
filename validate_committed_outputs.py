"""Validate committed aggregate outputs used in the manuscript."""
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parent
status = pd.read_csv(root / "data/processed/06_corrected_analysis_status.csv").iloc[0]
gdsc = pd.read_csv(root / "data/processed/06_gdsc_reproducibility.csv")

assert int(status["n_models"]) == 1139
assert int(status["n_expression_features"]) == 19130
assert int(status["n_testable_features"]) == 19123
assert int(status["n_fdr_significant"]) == 0
assert abs(float(status["minimum_fdr"]) - 0.4981034803694583) < 1e-10
assert status["legacy_network_status"] == "invalid_gene_alignment_do_not_use"

ln = gdsc.loc[gdsc.metric == "LN_IC50"].iloc[0]
auc = gdsc.loc[gdsc.metric == "AUC"].iloc[0]
assert int(ln.n_overlap) == 14
assert abs(float(ln.pearson_r) - 0.3258833074811336) < 1e-12
assert abs(float(ln.spearman_rho) - 0.42857142857142855) < 1e-12
assert abs(float(auc.pearson_r) - 0.2757980531796416) < 1e-12
print("Committed manuscript outputs validated.")
