# TCGA-CESC direct expression extraction. Dependency projection remains disabled
# until both model transferability and cross-platform calibration are established.
suppressPackageStartupMessages({library(TCGAbiolinks);library(SummarizedExperiment);library(tidyverse)})
dir.create("data/processed",recursive=TRUE,showWarnings=FALSE)
status<-read_csv("data/processed/02_model_status.csv",show_col_types=FALSE)
if(isTRUE(status$generalizable[[1]]))stop("A model was retained, but frozen TCGA calibration is not yet specified; projection halted.")
q<-GDCquery(project="TCGA-CESC",data.category="Transcriptome Profiling",data.type="Gene Expression Quantification",workflow.type="STAR - Counts")
GDCdownload(q);se<-GDCprepare(q);rn<-as.data.frame(rowData(se));x<-assay(se,"tpm_unstrand")
symbols<-as.character(rn$gene_name);extract_gene<-function(g){i<-which(toupper(symbols)==g);if(!length(i))return(rep(NA_real_,ncol(x)));apply(log2(x[i,,drop=FALSE]+1),2,median,na.rm=TRUE)}
out<-tibble(barcode=colnames(x),patient_id=substr(barcode,1,12),PRKDC=extract_gene("PRKDC"),BARD1=extract_gene("BARD1"))
hpv_file<-"data/raw/TCGA_CESC_validated_HPV_calls.tsv"
if(file.exists(hpv_file)){hpv<-read_tsv(hpv_file,show_col_types=FALSE);if(!all(c("patient_id","hpv_status")%in%names(hpv)))stop("HPV table must contain patient_id and hpv_status.");out<-left_join(out,hpv,by="patient_id")}else warning("Validated HPV calls absent; HPV-stratified analysis skipped.")
write_csv(out,"data/processed/03_tcga_cesc_direct_expression.csv")
