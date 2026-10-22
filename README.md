# PRKDC-cervical-cancer-validation
# PRKDC CRISPR Dependency in Cervical Cancer

This repository contains the computational analysis pipeline and figure-generation scripts for the manuscript: Running title: DNA-PK Dependency in Cervical Cancer 
## Overview
The study derives an elastic-net dependency model from DepMap CRISPR data and projects it into TCGA-CESC, GSE56363, GSE3578, and GSE224327 to evaluate its potential as a clinical response biomarker.

## Pipeline Execution
The scripts are numbered in the exact order they should be executed to reproduce the findings:
*   `00` to `01`: Data downloading and harmonization.
*   `02`: Model development and cross-validation.
*   `03` to `06b`: Clinical cohort projection (TCGA and GEO).
*   `07`: Falsification and stress tests.
*   `08` to `15`: Figure generation.

## Environment Setup
This project uses `renv` for R package management. Run `renv::restore()` in the project directory to install the exact package versions used in this analysis.
