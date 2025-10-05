library("tidyverse")
library("glue")

# Random forest fitting with year test set --------------------------------

fit_rf_year <- function(rf_data, testing_year = "2024", params, rf_seed, target, mode = "classification", file = NULL) {
  set.seed(rf_seed)
  
  if (!is.null(file) && file.exists(file)) {
    return(readRDS(file))
  }
  
  # Retain location for splitting, drop it for modeling
  params_no_location <- setdiff(params, "location")
  
  model_data <- rf_data |>
    filter(!is.na(.data[[target]])) |>
    select(imp, sample_date, all_of(target), all_of(params))  # Assumes 'params' does not include location
  
  training_data <- model_data |>
    filter(year(sample_date) != testing_year)
  
  rf_spec <- 
    rand_forest(
      mode = mode,
      trees = 500,
      mtry = tune(),
      min_n = tune()
    ) |>
    set_engine("ranger",
               importance = "impurity",
               sample.fraction = c(0.5, 0.5),
               keep.inbag = TRUE
    )
  
  tuning_grid <- extract_parameter_set_dials(rf_spec) |>
    finalize(training_data |> select(-c(imp, location, all_of(target)))) |>
    grid_regular(levels = 5)
  
  rf_workflow <- workflow() |>
    add_model(rf_spec) |>
    add_formula(as.formula(paste(target, "~ .")))
  
  tune_random_forests <- function(data) {
    data <- data |> select(-location)  # remove location before tuning
    vfold_cv(data, v = 5, repeats = 2) |>
      tune_grid(
        object = rf_workflow,
        grid = tuning_grid
      )
  }
  
  tuning_results <- training_data |>
    group_by(imp) |>
    nest() |>
    ungroup() |>
    mutate(tuning_results = future_map(data, tune_random_forests, .options = furrr_options(seed = rf_seed)))
  
  summarized_results <- tuning_results |>
    mutate(metrics = map(tuning_results, collect_metrics)) |>
    unnest(metrics) |>
    select(-c(data, tuning_results)) |>
    mutate(hyperparams = paste(mtry, min_n, sep = "_")) |>
    select(-c(mtry, min_n)) |>
    group_by(hyperparams, .metric) |>
    summarize(avg = mean(mean), .groups = "drop") |>
    separate(hyperparams, into = c("mtry", "min_n"), sep = "_")
  
  best_metric <- if (mode == "classification") "accuracy" else "rmse"
  
  best_tuning <- summarized_results |>
    filter(.metric == best_metric) |>
    slice_max(avg) |>
    mutate(across(-.metric, as.numeric)) |>
    slice_min(mtry)
  
  final_wf <- rf_workflow |> finalize_workflow(best_tuning)
  
  final_fit_function <- function(data) {
    training_data <- data |>
      filter(year(sample_date) != testing_year) |>
      select(-sample_date, -location)
    testing_data <- data |>
      filter(year(sample_date) == testing_year) |>
      select(-sample_date, -location)
    
    final_wf |>
      last_fit(make_splits(training_data, assessment = testing_data))
  }
  
  model_fits <- model_data |>
    nest_by(imp) |>
    ungroup() |>
    mutate(final_model = future_map(data, final_fit_function, .options = furrr_options(seed = rf_seed)))
  
  model_results <- list(
    tuning_summary = summarized_results,
    model_fits = model_fits
  )
  
  saveRDS(model_results, file)
  return(model_results)
}


# Random forest fitting with location test set ----------------------------


fit_rf_location <- function(rf_data, testing_location, params, rf_seed, target, mode = "classification", file = NULL) {
  set.seed(rf_seed)
  
  if (!is.null(file) && file.exists(file)) {
    return(readRDS(file))
  }
  
  # Retain location for splitting, drop it for modeling
  params_no_location <- setdiff(params, "location")
  
  model_data <- rf_data |>
    filter(!is.na(.data[[target]])) |>
    select(imp, sample_date, location, all_of(target), all_of(params))  # keep location for filtering
  
  training_data <- model_data |>
    filter(location != testing_location)
  
  rf_spec <- 
    rand_forest(
      mode = mode,
      trees = 500,
      mtry = tune(),
      min_n = tune()
    ) |>
    set_engine("ranger",
               importance = "impurity",
               sample.fraction = c(0.5, 0.5),
               keep.inbag = TRUE
    )
  
  tuning_grid <- extract_parameter_set_dials(rf_spec) |>
    finalize(training_data |> select(all_of(params_no_location))) |>
    grid_regular(levels = 5)
  
  rf_workflow <- workflow() |>
    add_model(rf_spec) |>
    add_formula(as.formula(paste(target, "~ .")))
  
  tune_random_forests <- function(data) {
    data <- data |> select(-location)  # drop location before tuning
    vfold_cv(data, v = 5, repeats = 2) |>
      tune_grid(
        object = rf_workflow,
        grid = tuning_grid
      )
  }
  
  tuning_results <- training_data |>
    group_by(imp) |>
    nest() |>
    ungroup() |>
    mutate(tuning_results = future_map(data, tune_random_forests, .options = furrr_options(seed = rf_seed)))
  
  summarized_results <- tuning_results |>
    mutate(metrics = map(tuning_results, collect_metrics)) |>
    unnest(metrics) |>
    select(-c(data, tuning_results)) |>
    mutate(hyperparams = paste(mtry, min_n, sep = "_")) |>
    select(-c(mtry, min_n)) |>
    group_by(hyperparams, .metric) |>
    summarize(avg = mean(mean), .groups = "drop") |>
    separate(hyperparams, into = c("mtry", "min_n"), sep = "_")
  
  if (nrow(summarized_results) == 0) {
    stop("Tuning failed: No results returned. Check input data and predictors.")
  }
  
  best_metric <- if (mode == "classification") "accuracy" else "rmse"
  
  best_tuning <- summarized_results |>
    filter(.metric == best_metric) |>
    slice_max(avg) |>
    mutate(across(-.metric, as.numeric)) |>
    slice_min(mtry)
  
  final_wf <- rf_workflow |> finalize_workflow(best_tuning)
  
  final_fit_function <- function(data) {
    training_data <- data |>
      filter(location != testing_location) |>
      select(-sample_date, -location)
    testing_data <- data |>
      filter(location == testing_location) |>
      select(-sample_date, -location)
    
    final_wf |>
      last_fit(make_splits(training_data, assessment = testing_data))
  }
  
  model_fits <- model_data |>
    nest_by(imp) |>
    ungroup() |>
    mutate(final_model = future_map(data, final_fit_function, .options = furrr_options(seed = rf_seed)))
  
  model_results <- list(
    tuning_summary = summarized_results,
    model_fits = model_fits
  )
  
  saveRDS(model_results, file)
  return(model_results)
}


# Function to get the variable importance metrics from the model fits
extract_importance <- function(model_fit){
  model_fit |> 
    extract_fit_parsnip() |> 
    vip::vi()
}


# Functions to plot ROC curves --------------------------------------------

# Function to compute sensitivity and specificity across all model probabilities:
get_sens_spec_lookup <- function(data){
  
  total <- data |> 
    count(passive_mclr_binary) |> 
    ungroup() |> 
    distinct() |> 
    pivot_wider(names_from = passive_mclr_binary, values_from = n)
  
  data |>
    select(.pred_detected, passive_mclr_binary) |>
    # arrange by the predicted probability
    arrange(desc(.pred_detected)) |> 
    # for a given probability how many TP and FP are there 
    # (the number of times a true positives at higher predicted probabilities)
    mutate(
      is_detected = (passive_mclr_binary == "detected"),
      tp = cumsum(is_detected),
      fp = cumsum(!is_detected),
      sensitivity = tp / total$detected,
      fpr = fp / total$not_detected,
      specificity = 1 - fpr) |> 
    select(sensitivity, specificity, fpr)
}  

# Function to get the sensitivity at a given specificity:
get_sensitivity <- function(x, data){
  data |>
    filter(specificity - x >= 0) |>
    slice_max(sensitivity) |>
    mutate(specificity = x, fpr = 1 - x) |>
    distinct()
}



# Logistic regression helper functions ------------------------------------

extract_attributes <- function(x) {
  attributes(x) %>%
    set_names(janitor::make_clean_names(names(.))) %>%
    as_tibble() %>%
    slice(1)
}

# function to scale and keep attributes when split by location
scale_and_get_attrs <- function(x) {
  s <- scale(x)
  list(
    values = as.numeric(s),
    center = attr(s, "scaled:center"),
    scale = attr(s, "scaled:scale"))
}
