# Protocol-compliant PRKDC dependency model development
suppressPackageStartupMessages({library(tidyverse); library(tidymodels); library(glmnet)})
set.seed(5591); dir.create("models", recursive=TRUE, showWarnings=FALSE); dir.create("data/processed", recursive=TRUE, showWarnings=FALSE)
input <- "data/processed/01_harmonized_depmap_discovery.csv"
if (!file.exists(input)) stop("Missing harmonized DepMap table: ", input)
df <- read_csv(input, show_col_types=FALSE) %>% distinct(ModelID, .keep_all=TRUE)
required <- c("ModelID","OncotreeLineage","PRKDC_effect")
if (!all(required %in% names(df))) stop("Missing columns: ", paste(setdiff(required,names(df)),collapse=", "))
drop_pattern <- "(_cn$|_mut$|PRKDC_effect|ModelID|OncotreeLineage|PrimaryOrMetastasis)"
predictors <- names(df)[vapply(df,is.numeric,logical(1)) & !grepl(drop_pattern,names(df),ignore.case=TRUE)]
if (length(predictors)<2) stop("Fewer than two transcriptomic predictors found.")
dat <- df %>% select(all_of(required),all_of(predictors)) %>% drop_na(PRKDC_effect,OncotreeLineage) %>% mutate(OncotreeLineage=factor(OncotreeLineage))
if (n_distinct(dat$OncotreeLineage)<5) stop("At least five lineages are required.")
outer <- group_vfold_cv(dat,group=OncotreeLineage,v=min(10L,n_distinct(dat$OncotreeLineage)))
spec <- linear_reg(penalty=tune(),mixture=tune()) %>% set_engine("glmnet")
rec <- recipe(PRKDC_effect~.,data=dat %>% select(-ModelID,-OncotreeLineage)) %>% step_zv(all_predictors()) %>% step_impute_median(all_numeric_predictors()) %>% step_normalize(all_numeric_predictors())
wf <- workflow() %>% add_recipe(rec) %>% add_model(spec)
grid <- crossing(penalty=10^seq(-5,0,length.out=25),mixture=seq(0,1,by=.25))
fit_fold <- function(split,id){
  train <- analysis(split) %>% select(-ModelID,-OncotreeLineage); test <- assessment(split)
  tuned <- tune_grid(wf,resamples=vfold_cv(train,v=5,strata=PRKDC_effect),grid=grid,metrics=metric_set(rmse,rsq))
  best <- select_best(tuned,metric="rmse"); fitted <- fit(finalize_workflow(wf,best),train)
  bind_cols(test %>% select(ModelID,OncotreeLineage,PRKDC_effect),predict(fitted,test)) %>% mutate(outer_fold=id,penalty=best$penalty,mixture=best$mixture)
}
oof <- map2_dfr(outer$splits,outer$id,fit_fold); baseline <- mean(dat$PRKDC_effect)
rmse_model <- rmse_vec(oof$PRKDC_effect,oof$.pred); rmse_null <- rmse_vec(oof$PRKDC_effect,rep(baseline,nrow(oof)))
r2_oof <- 1-sum((oof$PRKDC_effect-oof$.pred)^2)/sum((oof$PRKDC_effect-mean(oof$PRKDC_effect))^2)
generalizable <- is.finite(r2_oof) && r2_oof>0 && rmse_model<rmse_null
cervix_n <- sum(str_detect(as.character(dat$OncotreeLineage), regex("cerv", ignore_case=TRUE)))
decision <- case_when(cervix_n==0 ~ "not_estimable_cervical_no_observed_prkdc_screens",
                      generalizable ~ "retain_transcriptomic_model", TRUE ~ "select_null_model")
status <- tibble(n=nrow(dat),p=length(predictors),outer_folds=nrow(outer),cervix_n=cervix_n,
                 oof_rmse=rmse_model,null_rmse=rmse_null,oof_r2=r2_oof,generalizable=generalizable,
                 clinical_transfer_eligible=generalizable && cervix_n>=5,decision=decision)
write_csv(oof,"data/processed/02_nested_oof_predictions.csv"); write_csv(status,"data/processed/02_model_status.csv")
if(generalizable){
  tuned <- tune_grid(wf,resamples=vfold_cv(dat %>% select(-ModelID,-OncotreeLineage),v=10,strata=PRKDC_effect),grid=grid,metrics=metric_set(rmse))
  best <- select_best(tuned,metric="rmse"); saveRDS(fit(finalize_workflow(wf,best),dat %>% select(-ModelID,-OncotreeLineage)),"models/locked_prkdc_dependency_model.rds")
  saveRDS(list(predictors=predictors,decision=status$decision,parameters=best),"models/model_manifest.rds")
} else {
  unlink(c("models/locked_prkdc_dependency_model.rds","models/locked_preprocessing_recipe.rds"))
  saveRDS(list(intercept=baseline,predictors=character(),decision=status$decision),"models/model_manifest.rds")
}
print(status)
