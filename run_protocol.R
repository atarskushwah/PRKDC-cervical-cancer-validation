args<-commandArgs(trailingOnly=TRUE);mode<-if(length(args))args[[1]]else"core"
run<-function(script){message("\n>>> ",script);status<-system2("Rscript",script);if(status!=0)stop(script," failed with status ",status)}
run("00_protocol_preflight.R")
if(mode%in%c("core","all")){run("01_depmap_harmonization.R");run("02_model_development.R")}
if(mode=="all"){
  status<-read.csv("data/processed/02_model_status.csv")
  if(!isTRUE(status$generalizable[[1]])){run("03_tcga_cesc_download.R");run("04_gse56363_response_validation.R");run("05_gse3578_paired_dynamics.R");run("05b_extract_manuscript_metrics.R")}
  else message("Clinical projection intentionally paused pending frozen platform calibration.")
}
