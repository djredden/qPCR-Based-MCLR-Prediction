# this script fits the hierarchical models that account for censoring and missing data that were used 
# to understand the relationship between the variables deemed important in the random forest work and 
# probability of MC-LR in the passive samplers or the concentration of MC-LR in the SPE grab samples

library(tidyverse)
library(assertr)
library(brms)
library(tidybayes)
library(marginaleffects)
library(modelr)
library(glue)
library(ggtext)
library(patchwork)
library(pROC)
library(yardstick)
library(viridisLite)
source("code/stan-functions.R")
source("code/functions.R")

options(mc.cores = parallel::detectCores())

grab_col <- "#51c016"
passive_col <- "#4a0c6b"
year_colors <- c("#E69F00","#56B4E9", "#009E73")

theme_set(theme_bw(base_size = 17))

# Functions ---------------------------------------------------------------

extract_thresholds <- function(data, group_vars, predictor_col, prob_col) {
  data |> 
    group_by(across(all_of(group_vars))) |> 
    arrange(.data[[predictor_col]]) |> 
    summarize(
      crosses_50 = any(.data[[prob_col]] >= 0.5),
      mcye_at_50 = if (crosses_50) {
        approx(x = .data[[prob_col]], y = .data[[predictor_col]], xout = 0.5)$y
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

# Function to ensure that specific weeks are scaled appropriately for making model predictions
scale_weeks <- function(weeks, scaling_vals) {
  (weeks - scaling_vals$s_week_of_year$scaled_center) /
    scaling_vals$s_week_of_year$scaled_scale
}

# Read data ---------------------------------------------------------------

data_raw <- read_csv("data-cleaned/data-complete.csv") 

# Prep data for modelling -------------------------------------------------

data_prepped <- data_raw |> 
  # verify no replicate measurements in a week
  group_by(location, week, year) |> 
  mutate(n = n()) |> 
  verify(n == 1) |> 
  select(-n) |> 
  ungroup() |> 
  mutate(
    # create censoring indicators
    cens_mcye = if_else(is.na(mcye_ndaf) | mcye_ndaf > lod_mcye_ndaf, "none", "left"),
    cens_total_cyano = if_else(is.na(total_cyano) | total_cyano > lod_total_cyano, "none", "left"),
    
    # replace the censored values with their censoring limit 
    mcye = if_else(mcye_ndaf < lod_mcye_ndaf, lod_mcye_ndaf, mcye_ndaf),
    total_cyano = if_else(total_cyano < lod_total_cyano, lod_total_cyano, total_cyano),
    cens_ratio = if_else(cens_mcye == "left", "left", "none"),
    cens_ratio = if_else(is.na(mcye) | is.na(total_cyano), "none", cens_ratio),
   
     # log transform qPCR predictors
    mcye = log10(mcye),
    total_cyano = log10(total_cyano),
    
    # ratio between the two
    ratio = mcye / total_cyano,
    
    # created lagged versions of mcye and total cyano (and their cens and lod's)
    lag_mcye = lag(mcye),
    lag_total_cyano = lag(total_cyano),
    lag_ratio = lag(ratio),
    cens_lag_mcye = lag(cens_mcye),
    cens_lag_mcye = replace_na(cens_lag_mcye, "none"),
    cens_lag_total_cyano = lag(cens_total_cyano),
    cens_lag_total_cyano = replace_na(cens_lag_total_cyano, "none"),
    
    # create binary indicator for MC-LR on passive samplers
    passive_mclr_binary = if_else(passive_mclr > passive_mclr_lod, 1, 0),
    location = factor(location),
    
    # add an index term for including week in the model and allowing it to span the three years
    isoyear = isoyear(sample_date),  # use ISO year to avoid week overlaps
    week_of_year = isoweek(sample_date),
    .after = sample_date
  ) 

# variables to include in the model (based on RF results)
these_params <- c("week_of_year", "mcye", "lag_mcye", "tds_field", "cond_field", "chloride", 
                  "uv_254", "colour", "precip_7d_sum")

# Grand mean scaling of the data ------------------------------------------

# Grand-mean scaling
data_scaled <- data_prepped |> 
  mutate(
    across(
      .col = all_of(these_params),
      .fns = ~scale(.x)[,1],
      .names = "s_{.col}"),
    across(
      .col = all_of(these_params),
      .fns = ~scale(.x),
      .names = "attr_{.col}")) |> 
  filter(!is.na(passive_mclr_binary)) |> 
  as_tibble()

# Save scaling vals 
scaling_vals <- data_scaled |> 
  select(matches("^attr_")) |> 
  summarize(across(everything(), ~extract_attributes(.))) |> 
  pivot_longer(everything()) |> 
  mutate(name = str_replace(name, "attr", "s")) |> 
  unnest(value) |> 
  split(~name)

model_data <- data_scaled |> 
  select(-matches("^attr_")) |> 
  mutate(location_year = interaction(location, year), .after = location)

# Model fitting -----------------------------------------------------------

# Null model
null_bform <- bf(passive_mclr_binary ~ 1 + (1 | location) + (1 | year), family = "bernoulli")

# Without random slopes
no_slopes_bform <- 
  bf(passive_mclr_binary ~ s(s_week_of_year, bs = "fs", by = location_year) + 
       mi(s_mcye) + year +  (1 | location), family = "bernoulli") +
  bf(s_mcye | mi() ~ year + (1 | location), family = "student") +
  set_rescor(FALSE)

# With random slopes
mcye_loc_year_bform <- 
  bf(passive_mclr_binary ~ s(s_week_of_year, bs = "fs", by = location_year) +
       mi(s_mcye) +  year +  (1 + mi(s_mcye) | location), family = "bernoulli") + 
  bf(s_mcye | mi() ~ year + (1 | location), family = "student") +
  set_rescor(FALSE)

# With random slopes but without year-varying smooths
mcye_bform <- 
  bf(passive_mclr_binary ~ s(s_week_of_year, bs = "fs", by = location) +
       mi(s_mcye) +  year +  (1 + mi(s_mcye) | location), family = "bernoulli") + 
  bf(s_mcye | mi() ~ year + (1 | location), family = "student") +
  set_rescor(FALSE)

# All shared predictors
bform_shared <- bf(
  passive_mclr_binary ~ s(s_week_of_year, bs = "fs", by = location_year) +
    mi(s_mcye) + 
    mi(s_lag_mcye) +
    mi(s_tds_field) +
    mi(s_cond_field) +
    mi(s_chloride) + 
    mi(s_uv_254) +
    mi(s_colour) +
    mi(s_precip_7d_sum) +
    (1 + mi(s_mcye) + mi(s_lag_mcye) + mi(s_tds_field) + mi(s_cond_field) + mi(s_chloride) +
        mi(s_uv_254) + mi(s_colour) + mi(s_precip_7d_sum) | location),
  family = "bernoulli") +
  bf(s_mcye | mi() ~ year + (1 | location), family = "student")  +
  bf(s_lag_mcye | mi() ~ year + (1 | location), family = "student") + 
  bf(s_tds_field | mi() ~ mi(s_cond_field) + year + (1 | location), family = "student") +
  bf(s_uv_254 | mi() ~ mi(s_colour) + year + (1 | location), family = "student") +
  bf(s_colour | mi() ~ year + (1 | location), family = "student") +
  bf(s_chloride | mi() ~ mi(s_cond_field) + year + (1 | location), family = "student") +
  bf(s_cond_field | mi() ~ year + (1 | location), family = "student") +
  bf(s_precip_7d_sum | mi() ~  year + (1 | location), family = "student") +
  set_rescor(FALSE)

# Top passive predictors
bform_passive <- bf(
  passive_mclr_binary ~ s(s_week_of_year, bs = "fs", by = location_year) +
  mi(s_mcye) + 
    mi(s_lag_mcye) +
    mi(s_tds_field) +
    mi(s_cond_field) +
    mi(s_chloride) + 
    (1 + mi(s_mcye) + mi(s_lag_mcye) + mi(s_tds_field) + mi(s_cond_field) + mi(s_chloride) | location),
  family = "bernoulli") +
  bf(s_mcye | mi() ~ year + (1 | location), family = "student")  +
  bf(s_lag_mcye | mi() ~ year + (1 | location), family = "student") + 
  bf(s_tds_field | mi() ~ mi(s_cond_field) + year + (1 | location), family = "student") +
  bf(s_chloride | mi() ~ mi(s_cond_field) + year + (1 | location), family = "student") +
  bf(s_cond_field | mi() ~  year + (1 | location), family = "student") +
  set_rescor(FALSE)

# Priors
bprior <-
  set_prior("student_t(5, 0, 2.5)", class = "b", resp = "passivemclrbinary") +
  set_prior("normal(0, 1)", class = "Intercept", resp = "passivemclrbinary") 

# Fit the models
null_model <- fit_stan_model(
  file = "model-fits/null-model",
  seed = 1543,
  bform = null_bform,
  model_data = model_data,
  lower_bound = NULL,
  # overwrite = TRUE,
  adapt_delta = 0.99,
  max_treedepth = 12)

no_slopes_model <- fit_stan_model(
  file = "model-fits/no-slopes-model",
  seed = 1543,
  bform = no_slopes_bform,
  model_data = model_data,
  bpriors = bprior,
  var_xcens = c("s_mcye"),
  cens_ind = c("cens_mcye"),
  lower_bound = NULL,
  # overwrite = TRUE,
  adapt_delta = 0.99,
  max_treedepth = 12)

mcye_loc_year <- fit_stan_model(
  file = "model-fits/mcye-loc-year-model",
  seed = 1543,
  bform = mcye_loc_year_bform,
  model_data = model_data,
  bpriors = bprior,
  var_xcens = c("s_mcye"),
  cens_ind = c("cens_mcye"),
  lower_bound = NULL,
  # overwrite = TRUE,
  adapt_delta = 0.99,
  max_treedepth = 12)

mcye_model <- fit_stan_model(
  file = "model-fits/mcye-model",
  seed = 1543,
  bform = mcye_bform,
  model_data = model_data,
  bpriors = bprior,
  var_xcens = c("s_mcye"),
  cens_ind = c("cens_mcye"),
  lower_bound = NULL,
  # overwrite = TRUE,
  adapt_delta = 0.99,
  max_treedepth = 12)

# all shared predictors
shared_model <- fit_stan_model(
  file = "model-fits/shared-model",
  seed = 1543,
  bform = bform_shared,
  model_data = model_data,
  bpriors = bprior,
  var_xcens = c("s_mcye", "s_lag_mcye"),
  cens_ind = c("cens_mcye","cens_lag_mcye"),
  lower_bound = NULL,
  # overwrite = TRUE,
  adapt_delta = 0.999,
  max_treedepth = 15)

# passive top predictors
passive_model <- fit_stan_model(
  file = "model-fits/passive-model",
  seed = 1543,
  bform = bform_passive,
  model_data = model_data,
  bpriors = bprior,
  var_xcens = c("s_mcye", "s_lag_mcye"),
  cens_ind = c("cens_mcye","cens_lag_mcye"),
  lower_bound = NULL,
  # overwrite = TRUE,
  adapt_delta = 0.999,
  max_treedepth = 15)

# Model comparisons -------------------------------------------------------

# loo-cv (moment matching was needed for several models and was added to all for fair comparison)
# n.b. each model will be recompiled so this takes several minutes for all five models
loo_null  <- loo(null_model)
loo_no_slopes <- loo(no_slopes_model, resp = "passivemclrbinary")
loo_mcye  <- loo(mcye_model, resp = "passivemclrbinary")
loo_mcye_loc_year <- loo(mcye_loc_year, resp = "passivemclrbinary")
loo_shared <- loo(shared_model, resp = "passivemclrbinary")
loo_passive <- loo(passive_model, resp = "passivemclrbinary")

# compare the models
loo_list <- list(
  no_slopes_model = loo_no_slopes,
  mcye_shared     = loo_shared,
  mcye_loc_year   = loo_mcye_loc_year,
  passive_model   = loo_passive,
  null_model      = loo_null
)

loo_results <- loo_compare(loo_list)

model_comparison <- as_tibble(loo_results, rownames = "model") |>
  mutate(delta_looic = looic - min(looic))

# get the pareto-k warnings 
count_khat_exceedances <- function(loo_list, threshold = 0.7) {
  purrr::map_dfr(loo_list, ~{
    tibble(n_exceed = sum(.x$diagnostics$pareto_k > threshold, na.rm = TRUE))
  }, .id = "model")
}

# Get counts for k > 0.7 and k > 1
pareto_07 <- count_khat_exceedances(loo_list, threshold = 0.7) |>
  rename(n_khat_gt_0.7 = n_exceed)

pareto_1 <- count_khat_exceedances(loo_list, threshold = 1) |>
  rename(n_khat_gt_1 = n_exceed)

# Combine and join to model comparison table
pareto_summary <- pareto_07 |>
  left_join(pareto_1, by = "model")

model_comparison_augmented <- model_comparison |>
  left_join(pareto_summary, by = "model") |> 
  mutate(across(elpd_diff:delta_looic, ~ round(.x, digits = 2)))

write_csv(model_comparison_augmented, "tables/model_comparison_with_khat.csv")

# model checks ------------------------------------------------------------

# posterior predictive checks
yrep <- posterior_predict(mcye_loc_year, resp = "passivemclrbinary")
yrep <- as.matrix(yrep)

# mean detection rate per lake
ppc_mean <- bayesplot::ppc_stat_grouped(
  y = model_data$passive_mclr_binary,
  yrep = yrep,
  group = model_data$location,
  facet_args = list(scales = "free_y"),
  stat = "mean") +
  scale_x_continuous(
    breaks = c(0, 0.2, 0.4, 0.6),
    labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(
    values = c("y" = "black"),
    labels = c("Observed")) +
  scale_fill_manual(
    values = c("yrep" = "skyblue"),
    labels = c("Predicted")) +
  labs(
    tag = "b)",
    x = "Proportion of MC-LR detections", 
    color = NULL,
    fill = NULL) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.margin = margin(b = -10, unit = "pt"))

# check calibration across probability ranges
fitted_probs <- posterior_epred(mcye_loc_year, resp = "passivemclrbinary") |> colMeans()

observed_vs_predicted_probs <- model_data |> 
  mutate(
    pred_prob = fitted_probs,
    # cut into probability bins
    bin = cut(pred_prob, breaks = seq(0, 1, 0.1))
  ) |> 
  group_by(bin) |> 
  summarise(
    observed = mean(passive_mclr_binary), # proportion of observed detections in each bin
    predicted = mean(pred_prob), # the average predicted probability in the bin
    .groups = "drop"
  ) |> 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(size = 3) +
  geom_abline(linetype = "dashed") +
  labs(
    tag = "a)",
    x = "Predicted probability", 
    y = "Observed proportion") 

model_performance <- free(observed_vs_predicted_probs, side = "tb") + ppc_mean + plot_layout(ncol = 2, widths = c(3, 7)) &
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/model-performance.png",
       model_performance,
       height = 8, width = 14,
       dpi = 600)


# Posterior checks of variance and binary outcomes ------------------------

# first rename locations for plotting
pp_check_data <- model_data |> 
  mutate(location = fct_recode(location, 
                               "TRR" = "Tower Road Reservoir",
                               "TCR" = "Turtle Creek Reservoir",
                               "SC" = "Shubenacadie Canal"))

# standard deviation per lake
ppc_sd <- bayesplot::ppc_stat_grouped(
  y = pp_check_data$passive_mclr_binary,
  yrep = yrep,
  group = pp_check_data$location,
  facet_args = list(scales = "free_y"),
  stat = "sd") +
  scale_x_continuous(
    breaks = c(0, 0.2, 0.4, 0.6),
    labels = scales::number_format(accuracy = 0.1),
  ) +
  scale_color_manual(
    values = c("y" = "black"),
    labels = c("Observed sd")
  ) +
  scale_fill_manual(
    values = c("yrep" = "skyblue"),
    labels = c("Predicted sd")
  ) +
  labs(
    tag = "b)",
    x = "Standard deviation", 
    color = NULL,
    fill = NULL) +
  theme(legend.title = element_blank())


# binary outcome by lake

ppc_bars <- bayesplot::ppc_bars_grouped(
  y = pp_check_data$passive_mclr_binary,
  yrep = yrep,
  group = pp_check_data$location) +
  scale_x_continuous(
    breaks = c(0, 1),
    labels =c("Absent", "Present")
  ) + 
  scale_fill_manual(
    values = c("y" = "skyblue", "yrep" = "black"),
    labels = c("Observed count", "Predicted count")
  ) +
  scale_color_manual(
    values = c("yrep" = "black"),
    labels = c("Predicted count")
  ) +
  labs(
    tag = "a)",
    x = "MC-LR", 
    color = NULL,
    fill = NULL) +
  theme(
    legend.title = element_blank(),
    axis.ticks.x = element_blank())

model_checks_plot <- free(ppc_bars) + free(ppc_sd) & theme(legend.position = "top", plot.tag = element_text(face = "bold"))

ggsave("figures/si-figures/model-checks.png",
       model_checks_plot,
       height = 8, width = 14,
       dpi = 600)

# check ACFs
r <- resid(mcye_loc_year, resp = "passivemclrbinary")[,"Estimate"]
acf(r)

# Generate "thresholds" all locations and years ---------------------------

# Extract valid location–year combinations directly from training data
# This is required because location-year combos that don't exist (e.g. Shubie Canal in 2022) will cause errors
valid_df <- model_data |>
  distinct(location, year, location_year) 

valid_levels <- levels(mcye_loc_year$data$location_year)

# Sequences of mcye and weeks to predict over
s_mcye_seq <- seq(min(model_data$s_mcye, na.rm = TRUE),
                  max(model_data$s_mcye, na.rm = TRUE),
                  length.out = 200)

s_week_seq <- seq(min(model_data$s_week_of_year, na.rm = TRUE),
                  max(model_data$s_week_of_year, na.rm = TRUE),
                  length.out = 40)

# Target probabilities
target_probs <- c(0.5, 0.75, 0.95)

# File paths to save/load
thresholds_file <- "results/all_threshold_draws.rds"
preds_file <- "results/all_preds.rds"

# Flags for missing outputs
missing_thresholds <- !file.exists(thresholds_file)
missing_preds <- !file.exists(preds_file)

# Adding predicted draws across the grid defined by this combination of mcyE, week, year, and location
# is too big to store in memory so the predicted draws are added in chunks by looping over location and year
# loc_year is redefined within the loop because I was getting errors about new factor levels from predicted_draws. This
# approach ensures that only valid levels are included.

# N.B. This takes ~20 minutes to run

# Initialize containers
all_preds <- list()
all_threshold_draws <- list()

if (!missing_thresholds && !missing_preds) {
  message("Loading existing thresholds and predictions...")
  all_threshold_draws <- readRDS(thresholds_file)
  all_preds <- readRDS(preds_file)
  
} else {
  message("Generating missing predictions and/or thresholds...")
  
  for (i in seq_len(nrow(valid_df))) {
    loc <- valid_df$location[i]
    yr <- valid_df$year[i]
    loc_year <- as.character(valid_df$location_year[i])
    
    # Create prediction grid for this location × year
    grid_chunk <- expand_grid(
      location = loc,
      year = yr,
      s_week_of_year = s_week_seq,
      s_mcye = s_mcye_seq
    ) |>
      left_join(valid_df, by = c("location", "year")) |>
      mutate(
        location_year = factor(location_year, levels = valid_levels),
        mcye = s_mcye * scaling_vals$s_mcye$scaled_scale + scaling_vals$s_mcye$scaled_center,
        week_of_year = s_week_of_year * scaling_vals$s_week_of_year$scaled_scale + scaling_vals$s_week_of_year$scaled_center
      )
    
    # Posterior predictions
    epreds <- posterior_epred(
      mcye_loc_year,
      newdata = grid_chunk,
      resp = "passivemclrbinary",
      re_formula = NULL,
      ndraws = 4000,
      seed = 62
    )
    
    grid_chunk$mean_prob_draw <- colMeans(epreds)
    all_preds[[loc_year]] <- grid_chunk
    
    # Reshape to [draw × mcye × week]
    n_draws <- nrow(epreds) # 4000
    n_mcye  <- length(s_mcye_seq)
    n_week  <- length(s_week_seq)
    
    epreds_array <- array(epreds, dim = c(n_draws, n_mcye, n_week))
    mcye_vals    <- unique(grid_chunk$mcye)
    
    # Function to extract the thresholds
    extract_thresholds <- function(target_probs) {
      out_list <- vector("list", length(target_probs))
      
      for (pi in seq_along(target_probs)) {
        prob <- target_probs[pi]
        thr_mat <- matrix(NA_real_, nrow = n_draws, ncol = n_week)
        
        for (w in 1:n_week) {
          
          for (d in 1:n_draws) {
            p_vec <- epreds_array[d, , w]
            # interpolate mcyE values at target prob 
            # (see https://www.r-bloggers.com/2023/08/mastering-data-approximation-with-rs-approx-function/#google_vignette)
            thr_mat[d, w] <- if (length(unique(p_vec)) > 1)
              approx(x = p_vec, y = mcye_vals, xout = prob, rule = 2)$y else NA_real_
          }
        }
        # melt once
        out_list[[pi]] <- tibble(
          draw = rep(seq_len(n_draws), times = n_week),
          week_idx = rep(seq_len(n_week), each  = n_draws),
          threshold = as.vector(thr_mat),
          probability = prob
        )
      }
      bind_rows(out_list)
    }
    
    threshold_draws <- extract_thresholds(target_probs) |>
      mutate(
        location = loc,
        year = yr,
        week_of_year = s_week_seq[week_idx] * scaling_vals$s_week_of_year$scaled_scale +
          scaling_vals$s_week_of_year$scaled_center
      )
    
    all_threshold_draws[[loc_year]] <- threshold_draws
  }
  
  # Save outputs
  if (missing_thresholds) saveRDS(all_threshold_draws, thresholds_file)
  if (missing_preds) saveRDS(all_preds, preds_file)
}

# Combine predictions
all_year_preds <- bind_rows(all_preds)

# Summarize threshold draws
threshold_draws_df <- bind_rows(all_threshold_draws)

threshold_summary <- threshold_draws_df |>
  group_by(location, year, week_of_year, probability) |>
  summarise(
    median = median(threshold, na.rm = TRUE),
    q25 = quantile(threshold, 0.25, na.rm = TRUE),
    q75 = quantile(threshold, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    date = as.Date("2023-01-01") + lubridate::duration(7 * week_of_year, "days"),
    date = as.Date(date)
  )


# plot yearly probability curves by location ------------------------------

loc_yr_plot <- threshold_summary |>
  filter(probability == 0.5) |>
  ggplot(aes(x = date, y = median, color = as.factor(year), fill = as.factor(year))) +
  facet_wrap(vars(location)) +
  geom_line(linewidth = 1.25) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA, show.legend = FALSE) +
  scale_y_continuous(labels = scales::label_math()) +
  scale_fill_manual(values = year_colors) +
  scale_color_manual(values = year_colors) + 
  scale_x_date(
    date_breaks = "2 months",   # fewer major ticks
    date_labels = "%b",
    minor_breaks = NULL,
    expand = c(0, 0)
  ) +
  labs(
    col = NULL,
    fill = NULL,
    x = NULL,
    y = "*mcy*E &nbsp; (GC mL<sup>-1</sup>)") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.x = element_markdown(angle = 45, vjust = 1.4, hjust = 1.2),
    axis.title.y = element_markdown(),
    legend.margin =  margin(b= -10, unit="pt"),
    plot.margin = margin(l = 5, r = 15),
    strip.text = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2)
    )

ggsave("figures/p0.5-plot.png", loc_yr_plot,
       dev = "png", dpi = 600,
       width = 16, height = 8)


loc_yr_plot_0.75 <- threshold_summary |>
  filter(probability == 0.75) |>
  ggplot(aes(x = date, y = median, color = as.factor(year), fill = as.factor(year))) +
  facet_wrap(vars(location)) +
  geom_line(linewidth = 1.25) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA, show.legend = FALSE) +
  scale_y_continuous(labels = scales::label_math()) +
  scale_fill_manual(values = year_colors) +
  scale_color_manual(values = year_colors) +
  scale_x_date(
    date_breaks = "2 month",
    date_labels = "%b",
    expand = c(0, 0)) +
  labs(
    col = NULL,
    fill = NULL,
    x = NULL,
    y = "*mcy*E &nbsp; (GC mL<sup>-1</sup>)") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.x = element_markdown(angle = 45, vjust = 1.4, hjust = 1.2),
    axis.title.y = element_markdown(),
    legend.margin =  margin(b= -10, unit="pt"),
    plot.margin = margin(l = 5, r = 15),
    strip.text = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2)
  )


ggsave("figures/si-figures/p0.75-plot.png", loc_yr_plot_0.75,
       dev = "png", dpi = 600,
       width = 16, height = 8)

loc_yr_plot_0.95 <- threshold_summary |>
  filter(probability == 0.95) |>
  ggplot(aes(x = date, y = median, color = as.factor(year), fill = as.factor(year))) +
  facet_wrap(vars(location)) +
  geom_line(linewidth = 1.25) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA, show.legend = FALSE) +
  scale_y_continuous(labels = scales::label_math()) +
  scale_fill_manual(values = year_colors) +
  scale_color_manual(values = year_colors) +
  scale_x_date(
    date_breaks = "2 month",
    date_labels = "%b",
    expand = c(0, 0)) +
  labs(
    col = NULL,
    fill = NULL,
    x = NULL,
    y = "*mcy*E &nbsp; (GC mL<sup>-1</sup>)") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.x = element_markdown(angle = 45, vjust = 1.4, hjust = 1.2),
    axis.title.y = element_markdown(),
    legend.margin =  margin(b= -10, unit="pt"),
    plot.margin = margin(l = 5, r = 15),
    strip.text = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2)
  )


ggsave("figures/si-figures/p0.95-plot.png", loc_yr_plot_0.95,
       dev = "png", dpi = 600,
       width = 16, height = 8)


# Generate probability plot for graphical abstract ------------------------

# Define the week to plot
target_week <- 40
s_target_week <- scale_weeks(target_week, scaling_vals)

# Define the range of mcyE to plot over
mcye_range <- model_data |>
  summarize(min_mcye = 0, max_mcye = max(mcye, na.rm = TRUE))

mcye_seq <- seq(mcye_range$min_mcye, mcye_range$max_mcye, length.out = 200)

# Scale mcye
s_mcye_seq <- scale(mcye_seq, center = scaling_vals$s_mcye$scaled_center,
                    scale = scaling_vals$s_mcye$scaled_scale)[, 1]

# Define all combinations of location, year, and mcye
prediction_grid <- expand_grid(
  location = unique(model_data$location),
  year = unique(model_data$year),
  mcye = mcye_seq
) |>
  mutate(
    s_mcye = scale(mcye, center = scaling_vals$s_mcye$scaled_center,
                   scale = scaling_vals$s_mcye$scaled_scale)[, 1],
    s_week_of_year = s_target_week,
    week_of_year = target_week,
    location_year = paste(location, year, sep = ".")
  )

# Predict from the model
fitted_vals <- fitted(
  mcye_model,
  newdata = prediction_grid,
  resp = "passivemclrbinary",
  re_formula = NULL,
  scale = "response"
)

# Combine with input grid
plot_data <- bind_cols(prediction_grid, as_tibble(fitted_vals))

# Plot: Facet by lake, color by year
plot_data |> 
  filter(location %in% c("Lake Fletcher", "INP", "Lake Charles 2")) |> 
  mutate(location = fct_recode(location,
                               "Lake A" = "INP",
                               "Lake B" = "Lake Charles 2",
                               "Lake C" = "Lake Fletcher")) |> 
  ggplot(aes(x = mcye, y = Estimate, color = factor(year))) +
  geom_line(linewidth = 0.6) +
  facet_wrap(vars(location), ncol = 1) +
  scale_color_manual(values = year_colors) +
  scale_x_continuous(labels = scales::math_format()) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1),
    limits = c(0,1)
  ) +
  labs(
    col = NULL,
    x = "*mcy*E (GC mL<sup>-1</sup>)",
    y = "Probability of MC-LR"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_markdown(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 9),
    strip.text.x = element_text(face = "bold", size = 10),
    legend.position = "top",
    legend.margin = margin(b = -15, unit = "pt"),
    legend.text = element_text(size = 9)
  )

ggsave("extra-figures/ga-image.svg", width = 2.5, height = 3, units = "in")

