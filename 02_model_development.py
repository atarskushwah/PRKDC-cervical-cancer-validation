"""Protocol-compliant nested, lineage-grouped PRKDC dependency analysis.

This Python implementation is used when the execution environment lacks R. It
uses the same frozen decision rule as 02_model_development.R and writes the
same core output tables.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

SEED = 5591
ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

def find_one(stem: str) -> Path:
    hits = sorted(RAW.glob(f"{stem}*.csv"))
    if not hits:
        raise FileNotFoundError(f"No input matching {stem}*.csv")
    return hits[0]

crispr_path = find_one("CRISPRGeneEffect")
expr_path = find_one("OmicsExpressionTPMLogp1HumanProteinCodingGenes")
model_path = find_one("Model")

crispr = pd.read_csv(crispr_path, usecols=lambda c: c == "" or c.startswith("Unnamed:") or c == "PRKDC (5591)")
crispr.columns = ["ModelID", "PRKDC_effect"]
models = pd.read_csv(model_path, usecols=["ModelID", "OncotreeLineage"])
expr = pd.read_csv(expr_path, low_memory=False)
if "ModelID" not in expr.columns:
    raise ValueError("Expression file lacks ModelID")

metadata = {"", "SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"}
features = [c for c in expr.columns if c not in metadata]
dat = (crispr.merge(models, on="ModelID", how="inner")
             .merge(expr[["ModelID"] + features], on="ModelID", how="inner")
             .dropna(subset=["PRKDC_effect", "OncotreeLineage"])
             .drop_duplicates("ModelID"))
X = dat[features].apply(pd.to_numeric, errors="coerce").astype("float32")
y = dat["PRKDC_effect"].to_numpy(dtype="float64")
groups = dat["OncotreeLineage"].astype(str).to_numpy()

outer = GroupKFold(n_splits=min(5, np.unique(groups).size))
alphas = np.logspace(-4, 0, 8)
param_grid = {"model__alpha": alphas, "model__l1_ratio": [0.25, 0.5, 0.75]}
rows = []

for fold, (train, test) in enumerate(outer.split(X, y, groups), start=1):
    pipe = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", ElasticNet(max_iter=10000, selection="cyclic", random_state=SEED)),
    ])
    inner_groups = groups[train]
    inner = GroupKFold(n_splits=min(3, np.unique(inner_groups).size))
    search = GridSearchCV(pipe, param_grid, scoring="neg_root_mean_squared_error",
                          cv=inner, n_jobs=-1, refit=True, error_score="raise")
    search.fit(X.iloc[train], y[train], groups=inner_groups)
    pred = search.predict(X.iloc[test])
    null_pred = np.repeat(y[train].mean(), len(test))
    for idx, p, n in zip(test, pred, null_pred):
        rows.append({"ModelID": dat.iloc[idx]["ModelID"], "OncotreeLineage": groups[idx],
                     "PRKDC_effect": y[idx], ".pred": p, ".null_pred": n,
                     "outer_fold": fold, "penalty": search.best_params_["model__alpha"],
                     "mixture": search.best_params_["model__l1_ratio"]})

oof = pd.DataFrame(rows)
rmse = mean_squared_error(oof.PRKDC_effect, oof[".pred"]) ** 0.5
null_rmse = mean_squared_error(oof.PRKDC_effect, oof[".null_pred"]) ** 0.5
r2 = r2_score(oof.PRKDC_effect, oof[".pred"])
mae = mean_absolute_error(oof.PRKDC_effect, oof[".pred"])

rng = np.random.default_rng(SEED)
boot = []
for _ in range(2000):
    i = rng.integers(0, len(oof), len(oof))
    yy = oof.PRKDC_effect.to_numpy()[i]
    pp = oof[".pred"].to_numpy()[i]
    nn = oof[".null_pred"].to_numpy()[i]
    boot.append((r2_score(yy, pp), mean_squared_error(yy, pp) ** .5 - mean_squared_error(yy, nn) ** .5))
boot = np.asarray(boot)
generalizable = bool(np.isfinite(r2) and r2 > 0 and rmse < null_rmse)
cervix_n = int((np.char.find(np.char.lower(groups.astype(str)), "cerv") >= 0).sum())
# Cervical transfer cannot be evaluated if the outcome matrix contains no
# cervical models, irrespective of pan-lineage prediction metrics.
clinical_transfer_eligible = bool(generalizable and cervix_n >= 5)
if cervix_n == 0:
    decision = "not_estimable_cervical_no_observed_prkdc_screens"
elif not generalizable:
    decision = "select_null_model"
else:
    decision = "retain_transcriptomic_model"
status = pd.DataFrame([{
    "n": len(dat), "p": len(features), "lineages": np.unique(groups).size,
    "crispr_total_models": len(crispr), "prkdc_observed_models": int(crispr.PRKDC_effect.notna().sum()),
    "prkdc_missing_fraction": float(crispr.PRKDC_effect.isna().mean()), "cervix_n": cervix_n,
    "outer_folds": outer.n_splits, "oof_rmse": rmse, "null_rmse": null_rmse,
    "oof_mae": mae, "oof_r2": r2,
    "oof_r2_ci_low": np.nanpercentile(boot[:, 0], 2.5),
    "oof_r2_ci_high": np.nanpercentile(boot[:, 0], 97.5),
    "rmse_difference": rmse-null_rmse,
    "rmse_difference_ci_low": np.nanpercentile(boot[:, 1], 2.5),
    "rmse_difference_ci_high": np.nanpercentile(boot[:, 1], 97.5),
    "generalizable": generalizable, "clinical_transfer_eligible": clinical_transfer_eligible,
    "decision": decision,
}])
oof.to_csv(OUT / "02_nested_oof_predictions.csv", index=False)
status.to_csv(OUT / "02_model_status.csv", index=False)
(OUT / "02_run_manifest.json").write_text(json.dumps({
    "seed": SEED, "crispr_file": crispr_path.name, "expression_file": expr_path.name,
    "model_file": model_path.name, "outer_cv": "GroupKFold by OncotreeLineage",
    "inner_cv": "GroupKFold by OncotreeLineage", "decision_rule": "OOF R2 > 0 and RMSE < null RMSE"
}, indent=2))
print(status.to_string(index=False))
