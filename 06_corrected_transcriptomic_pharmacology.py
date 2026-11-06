"""Corrected lineage-adjusted expression and cervical NU7441 analyses."""
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, t as tdist

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

def bh_adjust(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    restored = np.empty(len(vals))
    restored[order] = np.minimum(adjusted, 1)
    out[ok] = restored
    return out

frozen = Path(os.environ.get("FROZEN_DISCOVERY_PATH", OUT / "01_harmonized_depmap_discovery.csv"))
if not frozen.exists():
    raise FileNotFoundError(f"Missing frozen discovery table: {frozen}")
dat = pd.read_csv(frozen, low_memory=False)
features = [c for c in dat.columns if c.endswith("_exp") and pd.api.types.is_numeric_dtype(dat[c])]
groups = dat.OncotreeLineage.fillna("Missing").astype(str)
y = dat.PRKDC_effect.to_numpy(float)
yr = y - pd.Series(y).groupby(groups).transform("mean").to_numpy()
X = dat[features].to_numpy(np.float32)
for _, idx in groups.groupby(groups).groups.items():
    idx = np.asarray(list(idx))
    X[idx] -= np.nanmean(X[idx], axis=0)
num = np.nansum(X * yr[:, None], axis=0)
den = np.sqrt(np.nansum(X * X, axis=0) * np.nansum(yr * yr))
r = np.divide(num, den, out=np.full_like(num, np.nan), where=den != 0)
df = len(dat) - groups.nunique() - 1
tstat = r * np.sqrt(df / (1 - r * r))
p = 2 * tdist.sf(np.abs(tstat), df)
correlations = pd.DataFrame({
    "Gene": [c[:-4] for c in features],
    "Pearson_R_lineage_adjusted": r,
    "P_value": p,
})
correlations["FDR_BH"] = bh_adjust(correlations.P_value)
correlations = correlations.sort_values("P_value", na_position="last")
correlations.to_csv(OUT / "06_lineage_adjusted_expression_correlations.csv", index=False)

gdsc_paths = {
    "GDSC1": os.environ.get("GDSC1_PATH", "data/raw/GDSC1_fitted_dose_response_24Jul22.xlsx"),
    "GDSC2": os.environ.get("GDSC2_PATH", "data/raw/GDSC2_fitted_dose_response_24Jul22.xlsx"),
}
frames = []
for dataset, path in gdsc_paths.items():
    path = Path(path)
    if not path.exists():
        continue
    d = pd.read_excel(path)
    keep = d.TCGA_DESC.eq("CESC") & d.DRUG_NAME.astype(str).str.upper().eq("NU7441")
    x = d.loc[keep, ["SANGER_MODEL_ID", "CELL_LINE_NAME", "COSMIC_ID", "DRUG_ID",
                     "DRUG_NAME", "LN_IC50", "AUC", "RMSE", "Z_SCORE"]].copy()
    x.insert(0, "DATASET", dataset)
    frames.append(x)
if frames:
    nu = pd.concat(frames, ignore_index=True)
    nu.to_csv(OUT / "06_gdsc_nu7441_cervical.csv", index=False)
    wide = nu.pivot(index=["SANGER_MODEL_ID", "CELL_LINE_NAME"], columns="DATASET",
                    values=["LN_IC50", "AUC"]).dropna()
    summary = []
    for metric in ("LN_IC50", "AUC"):
        a, b = wide[(metric, "GDSC1")], wide[(metric, "GDSC2")]
        pr, pp = pearsonr(a, b)
        sr, sp = spearmanr(a, b)
        summary.append({"metric": metric, "n_overlap": len(a), "pearson_r": pr,
                        "pearson_p": pp, "spearman_rho": sr, "spearman_p": sp})
    pd.DataFrame(summary).to_csv(OUT / "06_gdsc_reproducibility.csv", index=False)

testable = int(np.isfinite(correlations.P_value).sum())
audit = pd.DataFrame([{
    "n_models": len(dat),
    "n_lineages": groups.nunique(),
    "n_expression_features": len(features),
    "n_testable_features": testable,
    "n_fdr_significant": int((correlations.FDR_BH < .05).sum()),
    "minimum_fdr": correlations.FDR_BH.min(),
    "legacy_network_status": "invalid_gene_alignment_do_not_use",
}])
audit.to_csv(OUT / "06_corrected_analysis_status.csv", index=False)
print(audit.to_string(index=False))
