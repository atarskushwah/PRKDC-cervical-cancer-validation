"""Frozen PRKDC transcriptomic validation with lineage-separated resampling.

The cervical lineage is never used for feature selection, tuning, scaling,
imputation, or fitting. Null predictions use each training partition's mean.
"""
from pathlib import Path
import json
import os
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
RAW, OUT = ROOT / "data" / "raw", ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)


def load_harmonized():
    supplied = os.environ.get("FROZEN_DISCOVERY_PATH")
    paths = ([Path(supplied)] if supplied else []) + [
        OUT / "01_harmonized_depmap_discovery.csv",
        RAW / "01_harmonized_depmap_discovery.csv",
    ]
    path = next((p for p in paths if p.exists()), None)
    if path is None:
        raise FileNotFoundError("Frozen harmonized discovery table was not found")
    dat = pd.read_csv(path, low_memory=False).drop_duplicates("ModelID")
    required = {"ModelID", "OncotreeLineage", "PRKDC_effect"}
    if not required.issubset(dat.columns):
        raise ValueError(f"Missing columns: {sorted(required.difference(dat.columns))}")
    features = [c for c in dat.columns if c.endswith("_exp") and
                pd.api.types.is_numeric_dtype(dat[c])]
    if len(features) < 100:
        raise ValueError("Fewer than 100 numeric *_exp predictors found")
    return dat.dropna(subset=["PRKDC_effect", "OncotreeLineage"]).copy(), path, features


def make_pipeline():
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", ElasticNet(max_iter=20000, random_state=SEED, tol=1e-4)),
    ])


dat, input_path, features = load_harmonized()
is_cervix = dat.OncotreeLineage.astype(str).str.contains("cerv", case=False, na=False)
discovery = dat.loc[~is_cervix].reset_index(drop=True)
cervical = dat.loc[is_cervix].reset_index(drop=True)
X = discovery[features].astype("float32")
y = discovery.PRKDC_effect.to_numpy(float)
groups = discovery.OncotreeLineage.astype(str).to_numpy()

outer = GroupKFold(n_splits=min(5, np.unique(groups).size))
grid = {"model__alpha": np.logspace(-4, 0, 8),
        "model__l1_ratio": [0.25, 0.5, 0.75]}
rows = []
for fold, (train, test) in enumerate(outer.split(X, y, groups), 1):
    inner_groups = groups[train]
    inner = GroupKFold(n_splits=min(3, np.unique(inner_groups).size))
    search = GridSearchCV(make_pipeline(), grid,
                          scoring="neg_root_mean_squared_error", cv=inner,
                          n_jobs=-1, refit=True, error_score="raise")
    search.fit(X.iloc[train], y[train], groups=inner_groups)
    candidate = search.predict(X.iloc[test])
    null = np.repeat(y[train].mean(), len(test))
    for idx, pred, null_pred in zip(test, candidate, null):
        rows.append({"ModelID": discovery.iloc[idx].ModelID,
                     "OncotreeLineage": groups[idx],
                     "PRKDC_effect": y[idx], ".candidate_pred": pred,
                     ".null_pred": null_pred, "outer_fold": fold,
                     "penalty": search.best_params_["model__alpha"],
                     "mixture": search.best_params_["model__l1_ratio"]})

oof = pd.DataFrame(rows)
rmse = mean_squared_error(oof.PRKDC_effect, oof[".candidate_pred"]) ** .5
null_rmse = mean_squared_error(oof.PRKDC_effect, oof[".null_pred"]) ** .5
r2 = r2_score(oof.PRKDC_effect, oof[".candidate_pred"])
generalizable = bool(np.isfinite(r2) and r2 > 0 and rmse < null_rmse)
oof[".selected_pred"] = oof[".candidate_pred"] if generalizable else oof[".null_pred"]

# Candidate predictions are retained for an auditable external falsification
# test. The selected prediction remains the null when the frozen rule fails.
final_cv = GroupKFold(n_splits=min(5, np.unique(groups).size))
final_search = GridSearchCV(make_pipeline(), grid,
                            scoring="neg_root_mean_squared_error", cv=final_cv,
                            n_jobs=-1, refit=True, error_score="raise")
final_search.fit(X, y, groups=groups)
cervical_out = cervical[["ModelID", "OncotreeLineage", "PRKDC_effect"]].copy()
if len(cervical_out):
    cervical_out[".candidate_pred"] = final_search.predict(cervical[features].astype("float32"))
    cervical_out[".null_pred"] = y.mean()
    cervical_out[".selected_pred"] = (cervical_out[".candidate_pred"] if generalizable
                                      else cervical_out[".null_pred"])
    cerv_r2 = r2_score(cervical_out.PRKDC_effect, cervical_out[".candidate_pred"])
    cerv_rmse = mean_squared_error(cervical_out.PRKDC_effect,
                                   cervical_out[".candidate_pred"]) ** .5
    cerv_null = mean_squared_error(cervical_out.PRKDC_effect,
                                   cervical_out[".null_pred"]) ** .5
else:
    cerv_r2 = cerv_rmse = cerv_null = np.nan

rng, boot = np.random.default_rng(SEED), []
for _ in range(2000):
    i = rng.integers(0, len(oof), len(oof))
    yy, pp, nn = (oof.PRKDC_effect.to_numpy()[i],
                  oof[".candidate_pred"].to_numpy()[i],
                  oof[".null_pred"].to_numpy()[i])
    boot.append((r2_score(yy, pp), mean_squared_error(yy, pp) ** .5 -
                 mean_squared_error(yy, nn) ** .5))
boot = np.asarray(boot)
status = pd.DataFrame([{
    "n_total": len(dat), "n_discovery": len(discovery),
    "n_cervical_holdout": len(cervical), "p_expression": len(features),
    "discovery_lineages": np.unique(groups).size, "outer_folds": outer.n_splits,
    "oof_candidate_rmse": rmse, "oof_null_rmse": null_rmse,
    "oof_candidate_r2": r2,
    "oof_candidate_mae": mean_absolute_error(oof.PRKDC_effect, oof[".candidate_pred"]),
    "oof_r2_ci_low": np.nanpercentile(boot[:, 0], 2.5),
    "oof_r2_ci_high": np.nanpercentile(boot[:, 0], 97.5),
    "rmse_difference": rmse - null_rmse,
    "rmse_difference_ci_low": np.nanpercentile(boot[:, 1], 2.5),
    "rmse_difference_ci_high": np.nanpercentile(boot[:, 1], 97.5),
    "generalizable": generalizable,
    "decision": "retain_transcriptomic_model" if generalizable else "select_null_model",
    "cervical_candidate_rmse": cerv_rmse,
    "cervical_null_rmse": cerv_null,
    "cervical_candidate_r2": cerv_r2,
}])
oof.to_csv(OUT / "02_nested_oof_predictions.csv", index=False)
cervical_out.to_csv(OUT / "02_cervical_holdout_predictions.csv", index=False)
status.to_csv(OUT / "02_model_status.csv", index=False)
(OUT / "02_run_manifest.json").write_text(json.dumps({
    "seed": SEED, "input_file": input_path.name,
    "eligible_predictors": "numeric columns ending _exp only",
    "cervical_policy": "excluded from all discovery fitting and tuning",
    "outer_cv": "GroupKFold by OncotreeLineage",
    "inner_cv": "GroupKFold by OncotreeLineage",
    "null_policy": "training-fold mean",
    "decision_rule": "OOF R2 > 0 and candidate RMSE < fold-specific null RMSE",
}, indent=2))
print(status.to_string(index=False))
