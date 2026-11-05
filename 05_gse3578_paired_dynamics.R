# Paired GSE3578 analysis after collapsing technical replicates by patient/timepoint.
suppressPackageStartupMessages({library(GEOquery);library(tidyverse);library(org.Hs.eg.db);library(AnnotationDbi)})
dir.create("data/processed",recursive=TRUE,showWarnings=FALSE)
status<-read_csv("data/processed/02_model_status.csv",show_col_types=FALSE);if(isTRUE(status$generalizable[[1]]))stop("Frozen cross-platform calibration required.")
gse<-getGEO("GSE3578",GSEMatrix=TRUE)[[1]];pd<-pData(gse);tt<-as.character(pd$title)
meta<-tibble(geo_accession=rownames(pd),patient_id=str_extract(tt,"^[0-9]+"),timepoint=case_when(str_detect(tt,"_p")~"pre",str_detect(tt,"_d")~"during",TRUE~NA_character_))
x<-exprs(gse);if(max(x,na.rm=TRUE)>50)x<-log2(x+1);fd<-fData(gse);acc<-sub("\\..*","",trimws(sub(",.*","",as.character(fd$GB_LIST))))
sm<-AnnotationDbi::mapIds(org.Hs.eg.db,keys=unique(acc[acc!=""]),column="SYMBOL",keytype="ACCNUM",multiVals="first");sy<-unname(sm[acc])
gv<-function(g){i<-which(toupper(sy)==g);if(!length(i))return(rep(NA_real_,ncol(x)));apply(x[i,,drop=FALSE],2,median,na.rm=TRUE)}
long<-meta%>%mutate(PRKDC=gv("PRKDC"),BARD1=gv("BARD1"))%>%drop_na(patient_id,timepoint)%>%group_by(patient_id,timepoint)%>%summarise(across(c(PRKDC,BARD1),median,na.rm=TRUE),.groups="drop")
paired<-long%>%pivot_wider(names_from=timepoint,values_from=c(PRKDC,BARD1))%>%drop_na()
tests<-map_dfr(c("PRKDC","BARD1"),function(g){a<-paired[[paste0(g,"_pre")]];b<-paired[[paste0(g,"_during")]];w<-wilcox.test(b,a,paired=TRUE,exact=FALSE);tibble(gene=g,patients=length(a),median_change=median(b-a),wilcoxon_p=w$p.value)})%>%mutate(fdr=p.adjust(wilcoxon_p,"BH"),analysis="exploratory_paired_expression")
write_csv(long,"data/processed/05_gse3578_collapsed_expression.csv");write_csv(tests,"data/processed/05_gse3578_paired_tests.csv");print(tests)
