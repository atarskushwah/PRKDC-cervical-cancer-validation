# Frozen v1.0 reanalysis decision record

## Primary question

Can baseline transcriptomics predict PRKDC CRISPR dependency outside the lineage context used for training? The primary estimand is out-of-fold prediction performance in completely held-out cancer lineages.

## Corrections to v1.0

| Previous implementation | Why invalid | Frozen correction |
|---|---|---|
| Penalty fixed at 0.01 to force non-zero coefficients | Selects the desired answer after observing model collapse | Tune penalty and mixture inside each outer fold; permit the null model |
| Ordinary fit described as nested, lineage-grouped CV | Performance was not estimated out of sample | Outer lineage-grouped folds with inner hyperparameter tuning |
| Missing genes in GEO cohorts replaced with zero | Zero is a biological value, not missingness, and platform scales differ | No projection unless a transferable model and a frozen calibration procedure both exist |
| Technical replicates treated as independent samples | Inflates sample size and understates uncertainty | Collapse GSE3578 replicates within patient and timepoint before paired testing |
| Weak purity correlation called proof of biology | Failure to detect confounding is not validation | Report purity association only as a sensitivity analysis |
| Numeric metadata admitted as transcriptomic predictors | Allows non-expression fields to enter the model | Restrict predictors to numeric columns ending `_exp` |
| Inner resampling ignored lineage | Hyperparameters can exploit lineage structure | Group both inner and outer resampling by `OncotreeLineage` |
| Global outcome mean used as the null comparator | Uses information from held-out observations | Calculate the null mean separately within every training fold |
| Cervical models included during discovery fitting | Prevents a genuinely external cervical test | Exclude all cervical models from fitting, tuning and preprocessing |
| APOBR reported at `r = 1.0` while PRKDC was absent from the network | Indicates target-column or gene-label misalignment | Withdraw the legacy network and recompute expression associations from named columns |

## Decision rule

Retain the transcriptomic model only when lineage-held-out R² is positive and its RMSE is lower than the intercept-only comparator. Otherwise, freeze the null model, stop clinical projection, and report direct PRKDC/BARD1 expression analyses as exploratory rather than surrogate dependency validation.

The frozen historical table contains observed PRKDC CRISPR outcomes for 18 cervical models. These models form an untouched external test set. They are not used for feature filtering, imputation, scaling, tuning or fitting. The separate 26Q1 audit contains no observed cervical PRKDC screens and is retained only as a release-comparability warning; its missing values are never merged with or substituted for the frozen historical outcomes.

## Frozen corrected results

- Harmonized dataset: 1,139 models, comprising 1,121 non-cervical discovery models and 18 cervical holdout models.
- Eligible predictors: 16,312 numeric expression fields ending `_exp`.
- Lineage-grouped nested elastic net: candidate RMSE 0.125550; fold-specific null RMSE 0.125623; out-of-fold R² -0.00261. The null model is selected.
- Untouched cervical test: candidate RMSE 0.128264; discovery-mean null RMSE 0.126330; R² -0.03912.
- Lineage-adjusted univariate screen: zero of 16,312 transcripts significant at BH FDR <0.05; minimum FDR 0.4248.
- GDSC NU7441 reproducibility across 14 shared cervical models is weak: ln(IC50) Pearson r=0.326 and Spearman rho=0.429; AUC Pearson r=0.276 and Spearman rho=0.191.

These results do not establish absence of DNA-PK biology. They show that baseline bulk transcriptomics does not provide a reproducible predictor of PRKDC CRISPR dependency or NU7441 monotherapy response in the available cervical models.

## Missing-data boundaries

GSE151666 outcomes must come from a de-identified author-provided mapping; HPV and stage do not reconstruct outcomes. Validated TCGA HPV calls are required for HPV-stratified claims. The 80-patient methylation set may be analysed only as an explicitly labelled subset until cohort completeness is verified. Release-matched DepMap copy-number and mutation matrices are sensitivity inputs and must not be substituted silently.

## Statistical outputs

- DepMap: outer-fold predictions, RMSE, intercept-only RMSE, and out-of-fold R².
- Cervical holdout: candidate and selected-model predictions, RMSE versus the discovery-mean null, and holdout R².
- Transcriptome-wide audit: lineage-adjusted Pearson correlations with BH FDR; no ranking is interpreted as a validated co-dependency network when FDR is null.
- GDSC1/GDSC2: cervical NU7441 ln(IC50), AUC and RMSE, plus Pearson and Spearman cross-screen agreement.
- GSE56363: patient-level PRKDC/BARD1 medians, Wilcoxon tests, AUC, and BH-adjusted exploratory p-values.
- GSE3578: patient/timepoint replicate-collapsed values, paired Wilcoxon tests, median within-patient change, and BH adjustment.
- TCGA purity: Spearman correlation labelled as a confounding sensitivity analysis.
