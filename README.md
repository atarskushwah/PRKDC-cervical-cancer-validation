# PRKDC dependency and NU7441 reproducibility in cervical cancer models

This repository contains the analysis code and aggregate outputs for the manuscript "Baseline Transcriptomics Does Not Predict PRKDC Dependency or NU7441 Sensitivity in Cervical Cancer Models."

## Study design

The frozen harmonised cohort contains 1,139 models. All 18 cervical models were excluded before preprocessing and model development. The remaining 1,121 non-cervical models were analysed with nested lineage-grouped elastic-net regression and a fold-specific training-mean comparator. Clinical projection was permitted only if the candidate achieved positive pooled out-of-fold R-squared and lower RMSE than the null model. The candidate failed this gate, so patient projection was not an eligible analysis.

The completed orthogonal analyses comprise:

- External evaluation in the untouched 18-model cervical cohort.
- Lineage-adjusted transcriptome-wide association testing with false-discovery control.
- Cross-screen comparison of NU7441 response in 14 cervical models shared by GDSC1 and GDSC2.

The legacy APOBR-led network is withdrawn because its target and gene labels were misaligned. The frozen historical cohort is kept separate from the DepMap Public 26Q1 coverage audit.

## Reproducibility

Run `Rscript run_protocol.R core` for the DepMap decision. The corrected expression and pharmacology analysis is implemented in `06_corrected_transcriptomic_pharmacology.py`. GitHub Actions checks the committed decision record and GDSC reproducibility results on every push and pull request.

Controlled source datasets are not redistributed. See `PROTOCOL_V1_REANALYSIS.md` for data boundaries, analysis gates and source requirements.
