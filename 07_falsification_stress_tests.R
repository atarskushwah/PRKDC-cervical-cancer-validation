# Purity sensitivity analysis; absence of association does not establish validity.
suppressPackageStartupMessages({library(tidyverse);library(TCGAbiolinks)})
f<-list.files("data/processed",pattern="tcga.*predictions\\.csv$",full.names=TRUE,ignore.case=TRUE);if(!length(f))stop("No TCGA prediction file; purity analysis not run.")
p<-read_csv(f[[1]],show_col_types=FALSE);bc<-intersect(names(p),c("Patient_ID","ModelID","barcode","geo_accession"))[1]
if(is.na(bc)||!"Predicted_Dependency_Score"%in%names(p))stop("Incomplete TCGA prediction schema.")
p<-p%>%mutate(sample_barcode=substr(.data[[bc]],1,15));data("Tumor.purity")
q<-Tumor.purity%>%filter(Cancer.type=="CESC")%>%transmute(sample_barcode=substr(Sample.ID,1,15),CPE=as.numeric(CPE))%>%drop_na();d<-inner_join(p,q,by="sample_barcode")
if(nrow(d)<10)stop("Fewer than 10 matched samples.");ct<-cor.test(d$Predicted_Dependency_Score,d$CPE,method="spearman",exact=FALSE)
out<-tibble(n=nrow(d),rho=unname(ct$estimate),p=ct$p.value,interpretation="Sensitivity analysis only: weak association does not validate the score; strong association indicates possible purity confounding.")
write_csv(out,"data/processed/07_purity_sensitivity.csv");print(out)
