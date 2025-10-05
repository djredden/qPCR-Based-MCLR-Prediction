
library(tidyverse)
library(tidymodels)
library(furrr)
library(future)
source("code/functions.R")

plan(multisession) # enable parallelization of model tuning for random forest
R_RANGER_NUM_THREADS <- parallel::detectCores()

theme_set(theme_bw(base_size = 17))

# Read the data -----------------------------------------------------------

miss_ranger_imputation <- readRDS("model-fits/miss_ranger_imputation.rds")

# Prep for model ----------------------------------------------------------

n_datasets <- length(miss_ranger_imputation)

# Compile resulting datasets
missranger_results <- miss_ranger_imputation |> 
  set_names(paste0("imp_", 1:n_datasets)) |> 
  list_rbind(names_to = "imp") |> 
  mutate(
    month = month(sample_date),
    imp = str_remove(imp, "imp_"),
    imp = factor(imp, levels = 1:n_datasets),
    passive_mclr_binary = if_else(passive_mclr > passive_mclr_lod, "detected", "not_detected"),
    passive_mclr_binary = factor(passive_mclr_binary, levels = c("detected", "not_detected")),
    # create nutrient ratio parameters (TN and TOC are in mg/L, TP is µg/L)
    tn_tp_ratio = (tn / 14) / ((p / 1000) / 31),
    tc_tn_ratio = (toc / 12) /  (tn / 14),
    tc_tp_ratio = (toc / 12) / ((p / 1000) / 31),
    lag_mcye = lag(mcye_ndaf),
    lag_tc = lag(total_cyano)
  ) |> 
  relocate(passive_mclr_binary, .after = passive_mclr) 

# Declare parameters to include in the models -----------------------------

# full model
params <- c("location", "temp_field", "do_field", "cond_field", "ph_field",
            "tds_field", "al", "p", "fe", "chloride", "nitrate", "sulfate", "uv_254",
            "turbidity", "colour", "doc", "toc", "tn", "mcye_ndaf", "total_cyano",
            "max_temp_c", "min_temp_c", "mean_temp_c", "temp_7d_mean",
            "temp_3d_mean", "temp_7d_avgmax", "temp_3d_avgmax", "total_precip_mm", 
            "precip_7d_sum", "precip_3d_sum", "tn_tp_ratio", 
            "lag_mcye", "lag_tc")

# fit the models ----------------------------------------------------------

# full set of predictors holding out years
full_test_2024 <- fit_rf_year(missranger_results, testing_year = "2024", 
                              target = "passive_mclr_binary",
                              params = params, rf_seed = 101,
                              file = "model-fits/passive_2024.rds")

# full set of predictors holding out locations

full_test_b2 <- fit_rf_location(missranger_results, params = params, 
                                target = "passive_mclr_binary",
                                testing_location = "Lake Banook 2", rf_seed = 101,
                                file = "model-fits/passive_b2.rds")

full_test_p <- fit_rf_location(missranger_results, params = params, 
                                target = "passive_mclr_binary",
                                testing_location = "Penhorn Lake", rf_seed = 101,
                                file = "model-fits/passive_p.rds")

full_test_inp <- fit_rf_location(missranger_results, params = params,
                                 target = "passive_mclr_binary",
                                 testing_location = "INP", rf_seed = 101,
                                 file = "model-fits/passive_inp.rds")

full_test_fletcher <- fit_rf_location(missranger_results, params = params, 
                                      testing_location = "Lake Fletcher", 
                                      target = "passive_mclr_binary", rf_seed = 101,
                                      file = "model-fits/passive_fletcher.rds")

full_test_cunard <- fit_rf_location(missranger_results, params = params, 
                                    target = "passive_mclr_binary",
                                    testing_location = "Cunard Pond", rf_seed = 101,
                                    file = "model-fits/passive_cunard.rds")

full_test_sc <- fit_rf_location(missranger_results, params = params, 
                                target = "passive_mclr_binary",
                                testing_location = "Shubenacadie Canal", rf_seed = 101,
                                file = "model-fits/passive_sc.rds")

full_test_kearney <- fit_rf_location(missranger_results, params = params, 
                                     target = "passive_mclr_binary",
                                     testing_location = "Kearney Lake", rf_seed = 101,
                                     file = "model-fits/passive_kearney.rds")

full_test_lc1 <- fit_rf_location(missranger_results, params = params, 
                                 target = "passive_mclr_binary",
                                 testing_location = "Lake Charles 1", rf_seed = 101,
                                 file = "model-fits/passive_lc1.rds")

full_test_lc2 <- fit_rf_location(missranger_results, params = params, 
                                 target = "passive_mclr_binary",
                                 testing_location = "Lake Charles 2", rf_seed = 101,
                                 file = "model-fits/passive_lc2.rds")

full_test_trr <- fit_rf_location(missranger_results, params = params, 
                                 target = "passive_mclr_binary",
                                 testing_location = "Tower Road Reservoir", rf_seed = 101,
                                 file = "model-fits/passive_trr.rds")

full_test_tcr <- fit_rf_location(missranger_results, params = params, 
                                 target = "passive_mclr_binary",
                                 testing_location = "Turtle Creek Reservoir", rf_seed = 101,
                                 file = "model-fits/passive_tcr.rds")

# Process model fits ------------------------------------------------------

# Search environment for the model objects
model_names <- ls(pattern = "[full]_test_.*")

# Compile all the model fits and get the metrics and predictions
model_results <- setNames(lapply(mget(model_names), \(u) u[[2]]), model_names) |> 
  bind_rows(.id = "model") |> 
  mutate(
    metrics = map(final_model, collect_metrics),
    importance = map(final_model, extract_importance),
    predictions = map(final_model, collect_predictions),
    auc = map(predictions, ~ yardstick::roc_auc(data = .x, passive_mclr_binary, .pred_detected)),
    roc = map(predictions, ~ yardstick::roc_curve(data = .x, truth =  passive_mclr_binary, .pred_detected))
  )

# ROC curves --------------------------------------------------------------

# get AUC:
auc <- model_results |> 
  select(model, auc) |> 
  unnest(auc) |> 
  group_by(model) |> 
  summarise(
    min = quantile(.estimate, probs = 0.025),
    auc = median(.estimate),
    max = quantile(.estimate, probs = 0.975)
  ) |> 
  mutate(across(-model, ~ round(.x, 3))) |> 
  arrange(desc(auc)) |> 
  mutate(location = recode(model,
                           "full_test_sc" = "SC",
                           "full_test_lc1" = "Lake Charles 1",
                           "full_test_cunard" = "Cunard Pond",
                           "full_test_lc2" = "Lake Charles 2",
                           "full_test_kearney" = "Kearney Lake",
                           "full_test_fletcher" = "Lake Fletcher",
                           "full_test_p" = "Penhorn Lake",
                           "full_test_inp" = "INP",
                           "full_test_trr" = "TRR",
                           "full_test_tcr" = "TCR",
                           "full_test_b2" = "Lake Banook 2"))

write_csv(auc, "results/passive-auc.csv")

# fit curves:
specificity <- seq(0, 1, 0.01) # grid of specificities to compute sensitivities

roc_data <- model_results |> 
  select(model, imp, curve_data = roc) 

roc_plot_data <- roc_data |> 
  select(model, imp, curve_data) |> 
  unnest(curve_data) |> 
  group_by(model, imp, specificity) |> 
  summarise(
    sensitivity = median(sensitivity), .groups = "drop"
  ) |> 
  group_by(model, specificity) |> 
  summarise(
    lower = quantile(sensitivity, probs = 0.25),
    upper = quantile(sensitivity, probs = 0.75),
    sensitivity = median(sensitivity), .groups = "drop"
  ) |> 
  mutate(model_type = if_else(str_detect(model, "_2024$"), "Tested on 2024", "Tested on Location")) |> 
  arrange(model, specificity) |> 
  # averaging across models results in non-monotonicity in a few instances
  # use the cumulative minumum to enforce monotonic relationship between spec and sens
  group_by(model) |> 
  mutate(
    sensitivity = cummin(sensitivity),
    lower = cummin(lower),
    upper = cummin(upper)
    ) |> 
  ungroup() 
       
write_csv(roc_plot_data, "results/passive-roc-data.csv")    

# variable importance -----------------------------------------------------

var_imp_data_full <- model_results |>
  select(imp, model, importance) |> 
  unnest(importance) |> 
  janitor::clean_names() |> 
  mutate(
    model_type = if_else(str_detect(model, "_2024$"), "Tested on 2024", "Tested on Location")
  ) |> 
  group_by(model_type, model, variable) |> 
  # first get median importance of each variable over all the iterations of each model
  summarise(importance = median(importance), .groups = "drop") |> 
  # get the importance relative to the max for each variable within each model
  group_by(model_type, model) |> 
  mutate(
    rel_importance_within_models = 100 * importance / max(importance)
  ) |> 
  ungroup() |> 
  # summarise this over the models  
  group_by(model_type, variable) |> 
  summarise(
    avg_importance = median(rel_importance_within_models),
    lower = quantile(rel_importance_within_models, 0.25),
    upper = quantile(rel_importance_within_models, 0.75),
    .groups = "drop"
  ) |> 
  arrange(desc(avg_importance)) 

# write the data 
write_csv(var_imp_data_full, "results/var-imp-passive.csv")
