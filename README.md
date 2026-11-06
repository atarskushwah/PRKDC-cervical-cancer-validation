# PRKDC-cervical-cancer-validation
# PRKDC CRISPR Dependency in Cervical Cancer

This repository contains the computational analysis pipeline and figure-generation scripts for the manuscript: Running title: DNA-PK Dependency in Cervical Cancer 
## Overview
The study tests whether a transcriptomic model of DepMap PRKDC CRISPR dependency generalises across held-out cancer lineages. Clinical projection is conditional on successful nested validation; model collapse is retained and reported as a valid negative result.

## Pipeline Execution
The scripts are numbered in the exact order they should be executed to reproduce the findings:
*   `00` to `01`: Data downloading and harmonization.
*   `00_protocol_preflight.R`: machine-readable input audit; missing required data stop the run.
*   `02`: nested lineage-grouped model development with an intercept-only comparator.
*   `04` and `05`: patient-level clinical analyses used when the DepMap model selects the null.
*   `06`: corrected lineage-adjusted transcriptomic correlations and GDSC cervical NU7441 reproducibility.
*   `07`: Falsification and stress tests.
*   `08` to `15`: Figure generation.

## Environment Setup
This project uses `renv` for R package management. Run `renv::restore()` in the project directory to install the exact package versions used in this analysis.

Run `Rscript run_protocol.R core` for the DepMap decision or `Rscript run_protocol.R all` for eligible clinical analyses. See [PROTOCOL_V1_REANALYSIS.md](PROTOCOL_V1_REANALYSIS.md) for the frozen decision rules and remaining data boundaries.

The frozen historical analysis contains 1,139 models and an untouched 18-model cervical test set. It must not be combined with the separate DepMap 26Q1 coverage audit, in which cervical PRKDC measurements are absent. The legacy APOBR-led network is withdrawn because its target/gene labels are misaligned; `06_corrected_transcriptomic_pharmacology.py` generates the replacement FDR-controlled expression audit.
