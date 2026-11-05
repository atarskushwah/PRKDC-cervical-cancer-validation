# Release-matched DepMap harmonisation with explicit ProfileID-to-ModelID mapping.
suppressPackageStartupMessages({library(tidyverse);library(data.table)})
dir.create("data/processed",recursive=TRUE,showWarnings=FALSE)
need<-c("Model.csv","OmicsProfiles.csv","CRISPRGeneEffect.csv","OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
missing<-need[!file.exists(file.path("data/raw",need))];if(length(missing))stop("Missing: ",paste(missing,collapse=", "))
read_dm<-function(file){as_tibble(fread(file.path("data/raw",file),data.table=FALSE,check.names=FALSE))}
models<-read_dm("Model.csv");profiles<-read_dm("OmicsProfiles.csv")
if("DepMap_ID"%in%names(models))models<-rename(models,ModelID=DepMap_ID)
if(!"ModelID"%in%names(models))stop("Model.csv lacks ModelID/DepMap_ID.")

map_ids<-function(tab,label){
  id<-names(tab)[1];names(tab)[1]<-"source_id"
  if(all(grepl("^ACH-",tab$source_id)))return(rename(tab,ModelID=source_id))
  pid<-intersect(names(profiles),c("ProfileID","ProfileId","profile_id"))[1]
  mid<-intersect(names(profiles),c("ModelID","DepMap_ID"))[1]
  if(is.na(pid)||is.na(mid))stop("OmicsProfiles.csv lacks ProfileID/ModelID mapping for ",label)
  map<-profiles%>%transmute(source_id=as.character(.data[[pid]]),ModelID=as.character(.data[[mid]]))%>%distinct()
  out<-left_join(tab%>%mutate(source_id=as.character(source_id)),map,by="source_id")
  if(mean(is.na(out$ModelID))>.05)stop(label," mapping failure exceeds 5%.")
  out%>%select(-source_id)%>%relocate(ModelID)
}
crispr<-map_ids(read_dm("CRISPRGeneEffect.csv"),"CRISPR");expr<-map_ids(read_dm("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),"expression")
prkdc_col<-function(x){hit<-grep("^PRKDC( |\\(|$)|5591",names(x),ignore.case=TRUE,value=TRUE);if(!length(hit))stop("PRKDC column not found.");hit[[1]]}
outcome<-crispr%>%transmute(ModelID,PRKDC_effect=as.numeric(.data[[prkdc_col(crispr)]]))%>%drop_na()
out<-outcome%>%left_join(models,by="ModelID")%>%left_join(expr,by="ModelID")

# Optional release-matched sensitivity covariates; absence is recorded by preflight.
cnfile<-"data/raw/OmicsCNGeneWGS.csv"
if(file.exists(cnfile)){cn<-map_ids(read_dm(basename(cnfile)),"copy number");out<-out%>%left_join(cn%>%transmute(ModelID,PRKDC_cn=as.numeric(.data[[prkdc_col(cn)]])),by="ModelID")}
for(kind in c("Hotspot","Damaging")){
  f<-paste0("OmicsSomaticMutationsMatrix",kind,".csv");if(file.exists(file.path("data/raw",f))){m<-map_ids(read_dm(f),tolower(kind));col<-prkdc_col(m);v<-m%>%transmute(ModelID,value=as.numeric(.data[[col]]));names(v)[2]<-paste0("PRKDC_mut_",tolower(kind));out<-left_join(out,v,by="ModelID")}
}
if(anyDuplicated(out$ModelID))stop("Harmonisation produced duplicate ModelIDs.")
write_csv(out,"data/processed/01_harmonized_depmap_discovery.csv")
write_csv(tibble(models=nrow(out),expression_features=ncol(expr)-1),"data/processed/01_harmonization_summary.csv")
