# Frozen PRKDC transcriptomic validation with a completely held-out cervix lineage.
suppressPackageStartupMessages({library(tidyverse); library(tidymodels); library(glmnet)})
set.seed(5591)
dir.create("models", recursive=TRUE, showWarnings=FALSE)
dir.create("data/processed", recursive=TRUE, showWarnings=FALSE)
input <- Sys.getenv("FROZEN_DISCOVERY_PATH", "data/processed/01_harmonized_depmap_discovery.csv")
if (!file.exists(input)) stop("Missing frozen harmonized table: ", input)
df <- read_csv(input, show_col_types=FALSE) %>% distinct(ModelID, .keep_all=TRUE)
required <- c("ModelID", "OncotreeLineage", "PRKDC_effect")
if (!all(required %in% names(df))) stop("Missing columns: ", paste(setdiff(required,names(df)),collapse=", "))

# Numeric metadata are present in the frozen table. Only gene-expression
# columns with the explicit _exp suffix are eligible predictors.
predictors <- names(df)[endsWith(names(df), "_exp") & vapply(df,is.numeric,logical(1))]
if (length(predictors)<100) stop("Fewer than 100 numeric *_exp predictors found.")
dat <- df %>% select(all_of(required),all_of(predictors)) %>%
  drop_na(PRKDC_effect,OncotreeLineage) %>% mutate(OncotreeLineage=as.character(OncotreeLineage))
cervical <- dat %>% filter(str_detect(OncotreeLineage,regex("cerv",ignore_case=TRUE)))
discovery <- dat %>% filter(!str_detect(OncotreeLineage,regex("cerv",ignore_case=TRUE))) %>%
  mutate(OncotreeLineage=factor(OncotreeLineage))
if (n_distinct(discovery$OncotreeLineage)<5) stop("At least five non-cervical lineages are required.")

spec <- linear_reg(penalty=tune(),mixture=tune()) %>% set_engine("glmnet")
rec <- recipe(PRKDC_effect~.,data=discovery %>% select(PRKDC_effect,all_of(predictors))) %>%
  step_zv(all_predictors()) %>% step_impute_median(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors())
wf <- workflow() %>% add_recipe(rec) %>% add_model(spec)
grid <- crossing(penalty=10^seq(-4,0,length.out=8),mixture=c(.25,.5,.75))
outer <- group_vfold_cv(discovery,group=OncotreeLineage,v=min(5L,n_distinct(discovery$OncotreeLineage)))

fit_fold <- function(split,id){
  train_full <- analysis(split); test <- assessment(split)
  train_xy <- train_full %>% select(PRKDC_effect,all_of(predictors))
  inner <- group_vfold_cv(train_full,group=OncotreeLineage,v=min(3L,n_distinct(train_full$OncotreeLineage)))
  # Preserve grouped split indices while presenting only expression fields to the workflow.
  inner$splits <- map(inner$splits,function(s){
    make_splits(list(analysis=as.integer(s$in_id),assessment=as.integer(setdiff(seq_len(nrow(train_full)),s$in_id))),train_xy)
  })
  tuned <- tune_grid(wf,resamples=inner,grid=grid,metrics=metric_set(rmse))
  best <- select_best(tuned,metric="rmse")
  fitted <- fit(finalize_workflow(wf,best),train_xy)
  candidate <- predict(fitted,test %>% select(all_of(predictors)))$.pred
  null <- rep(mean(train_full$PRKDC_effect),nrow(test))
  test %>% select(ModelID,OncotreeLineage,PRKDC_effect) %>%
    mutate(.candidate_pred=candidate,.null_pred=null,outer_fold=id,
           penalty=best$penalty,mixture=best$mixture)
}
oof <- map2_dfr(outer$splits,outer$id,fit_fold)
rmse_model <- rmse_vec(oof$PRKDC_effect,oof$.candidate_pred)
rmse_null <- rmse_vec(oof$PRKDC_effect,oof$.null_pred)
r2_oof <- 1-sum((oof$PRKDC_effect-oof$.candidate_pred)^2)/sum((oof$PRKDC_effect-mean(oof$PRKDC_effect))^2)
generalizable <- is.finite(r2_oof) && r2_oof>0 && rmse_model<rmse_null
oof <- oof %>% mutate(.selected_pred=if (generalizable) .candidate_pred else .null_pred)

# Final candidate fit excludes every cervical model. Candidate predictions are
# retained for falsification even if the frozen rule selects the null.
final_inner <- group_vfold_cv(discovery,group=OncotreeLineage,v=min(5L,n_distinct(discovery$OncotreeLineage)))
final_xy <- discovery %>% select(PRKDC_effect,all_of(predictors))
final_inner$splits <- map(final_inner$splits,function(s){
  make_splits(list(analysis=as.integer(s$in_id),assessment=as.integer(setdiff(seq_len(nrow(discovery)),s$in_id))),final_xy)
})
final_tuned <- tune_grid(wf,resamples=final_inner,grid=grid,metrics=metric_set(rmse))
best <- select_best(final_tuned,metric="rmse")
final_fit <- fit(finalize_workflow(wf,best),final_xy)
cervical_out <- cervical %>% select(ModelID,OncotreeLineage,PRKDC_effect) %>%
  mutate(.candidate_pred=predict(final_fit,cervical %>% select(all_of(predictors)))$.pred,
         .null_pred=mean(discovery$PRKDC_effect),
         .selected_pred=if (generalizable) .candidate_pred else .null_pred)
cerv_r2 <- 1-sum((cervical_out$PRKDC_effect-cervical_out$.candidate_pred)^2)/sum((cervical_out$PRKDC_effect-mean(cervical_out$PRKDC_effect))^2)

status <- tibble(n_total=nrow(dat),n_discovery=nrow(discovery),n_cervical_holdout=nrow(cervical),
  p_expression=length(predictors),discovery_lineages=n_distinct(discovery$OncotreeLineage),
  outer_folds=nrow(outer),oof_candidate_rmse=rmse_model,oof_null_rmse=rmse_null,
  oof_candidate_r2=r2_oof,generalizable=generalizable,
  decision=if_else(generalizable,"retain_transcriptomic_model","select_null_model"),
  cervical_candidate_rmse=rmse_vec(cervical_out$PRKDC_effect,cervical_out$.candidate_pred),
  cervical_null_rmse=rmse_vec(cervical_out$PRKDC_effect,cervical_out$.null_pred),
  cervical_candidate_r2=cerv_r2)
write_csv(oof,"data/processed/02_nested_oof_predictions.csv")
write_csv(cervical_out,"data/processed/02_cervical_holdout_predictions.csv")
write_csv(status,"data/processed/02_model_status.csv")
saveRDS(list(intercept=mean(discovery$PRKDC_effect),predictors=character(),decision=status$decision),"models/model_manifest.rds")
print(status)
