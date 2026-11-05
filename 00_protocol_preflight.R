# Validate inputs without downloading, imputing, or substituting datasets.
suppressPackageStartupMessages(library(tidyverse))
dir.create("data/processed",recursive=TRUE,showWarnings=FALSE)
inputs<-tribble(
~stage,~file,~required,~purpose,
"depmap","data/raw/Model.csv",TRUE,"cell-line lineage metadata",
"depmap","data/raw/OmicsProfiles.csv",TRUE,"omics profile-to-model mapping",
"depmap","data/raw/CRISPRGeneEffect.csv",TRUE,"PRKDC CRISPR outcome",
"depmap","data/raw/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",TRUE,"transcriptomic predictors",
"copy_number","data/raw/OmicsCNGeneWGS.csv",FALSE,"copy-number sensitivity analysis",
"mutation","data/raw/OmicsSomaticMutationsMatrixHotspot.csv",FALSE,"release-matched hotspot covariates",
"mutation","data/raw/OmicsSomaticMutationsMatrixDamaging.csv",FALSE,"release-matched damaging covariates",
"gse56363","data/raw/GSE56363_series_matrix.txt.gz",FALSE,"response validation cache",
"gse3578","data/raw/GSE3578_series_matrix.txt.gz",FALSE,"paired dynamics cache",
"gse151666","data/raw/GSE151666_outcome_mapping.tsv",FALSE,"de-identified patient outcomes; never infer",
"tcga","data/raw/TCGA_CESC_validated_HPV_calls.tsv",FALSE,"validated HPV subgroup analysis",
"tcga","data/raw/TCGA_CESC_methylation_manifest.tsv",FALSE,"methylation cohort completeness audit")
report<-inputs%>%mutate(exists=file.exists(file),bytes=if_else(exists,as.numeric(file.info(file)$size),NA_real_),status=case_when(exists~"available",required~"BLOCKING",TRUE~"optional_or_stage_specific"))
write_csv(report,"data/processed/00_input_status.csv");print(report,n=Inf)
if(any(report$required&!report$exists))quit(status=2)
