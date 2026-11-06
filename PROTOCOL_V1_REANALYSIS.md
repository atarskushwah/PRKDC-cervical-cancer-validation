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

## Decision rule

Retain the transcriptomic model only when lineage-held-out R² is positive and its RMSE is lower than the intercept-only comparator. Otherwise, freeze the null model, stop clinical projection, and report direct PRKDC/BARD1 expression analyses as exploratory rather than surrogate dependency validation.

Clinical transfer additionally requires observed PRKDC CRISPR outcomes in cervical models. If the release contains no observed cervical PRKDC screens, cervical transferability is classified as not estimable rather than negative or positive.

## Missing-data boundaries

GSE151666 outcomes must come from a de-identified author-provided mapping; HPV and stage do not reconstruct outcomes. Validated TCGA HPV calls are required for HPV-stratified claims. The 80-patient methylation set may be analysed only as an explicitly labelled subset until cohort completeness is verified. Release-matched DepMap copy-number and mutation matrices are sensitivity inputs and must not be substituted silently.

## Statistical outputs

- DepMap: outer-fold predictions, RMSE, intercept-only RMSE, and out-of-fold R².
- GSE56363: patient-level PRKDC/BARD1 medians, Wilcoxon tests, AUC, and BH-adjusted exploratory p-values.
- GSE3578: patient/timepoint replicate-collapsed values, paired Wilcoxon tests, median within-patient change, and BH adjustment.
- TCGA purity: Spearman correlation labelled as a confounding sensitivity analysis.
