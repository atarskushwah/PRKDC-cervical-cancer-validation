# Direct patient-level GSE56363 analysis; null DepMap models are never projected.
suppressPackageStartupMessages({library(GEOquery);library(tidyverse);library(pROC)})
dir.create("data/processed",recursive=TRUE,showWarnings=FALSE)
status<-read_csv("data/processed/02_model_status.csv",show_col_types=FALSE)
if(isTRUE(status$generalizable[[1]])) stop("Frozen cross-platform calibration is required before projection.")
gse<-getGEO("GSE56363",GSEMatrix=TRUE)[[1]]; pd<-pData(gse)
txt<-apply(pd,1,paste,collapse=" "); response<-case_when(str_detect(tolower(txt),"non.?complete|_ncr_|partial|progressive")~0L,str_detect(tolower(txt),"complete response|_cr_")~1L,TRUE~NA_integer_)
fd<-fData(gse); sc<-grep("gene.?symbol|symbol",names(fd),ignore.case=TRUE,value=TRUE)[1]; if(is.na(sc)) stop("No gene-symbol annotation.")
x<-exprs(gse); if(max(x,na.rm=TRUE)>50)x<-log2(x+1); symbols<-sub("///.*|,.*","",as.character(fd[[sc]]))
gv<-function(g){i<-which(toupper(symbols)==g);if(!length(i))return(rep(NA_real_,ncol(x)));apply(x[i,,drop=FALSE],2,median,na.rm=TRUE)}
res<-tibble(geo_accession=colnames(x),response=response,PRKDC=gv("PRKDC"),BARD1=gv("BARD1"))%>%drop_na(response)
if(n_distinct(res$response)!=2)stop("Both response groups were not recovered.")
tests<-map_dfr(c("PRKDC","BARD1"),function(g){v<-res[[g]];w<-wilcox.test(v~res$response,exact=FALSE);ro<-roc(res$response,v,quiet=TRUE);tibble(gene=g,n=sum(is.finite(v)),median_complete=median(v[res$response==1],na.rm=TRUE),median_noncomplete=median(v[res$response==0],na.rm=TRUE),wilcoxon_p=w$p.value,auc=as.numeric(auc(ro)))})%>%mutate(fdr=p.adjust(wilcoxon_p,"BH"),analysis="exploratory_direct_expression")
write_csv(res,"data/processed/04_gse56363_patient_expression.csv");write_csv(tests,"data/processed/04_gse56363_direct_tests.csv");print(tests)
