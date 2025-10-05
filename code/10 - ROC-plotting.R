library(tidyverse)
library(glue)
library(ggtext)
library(patchwork)

grab_col <- "#51c016"
passive_col <- "#4a0c6b"

theme_set(theme_bw(base_size = 17))

# Functions ---------------------------------------------------------------

# Create the plots
var_imp_plotting <- function(df){
  
  df |> 
    ggplot(aes(x = avg_importance, y = variable, alpha = shared)) +
    geom_point(size = 3, col = "grey15") +
    geom_errorbar(aes(xmin = lower, xmax = upper, alpha = shared),
                  width = 0.2,  col = "grey15", linewidth = 1.2) +
    scale_alpha_manual(values = c(0.8, 0.2)) +
    scale_x_continuous(limits = c(30, 100)) +
    labs(
      x = "Relative importance (%)", 
      y = NULL) +
    theme(
      legend.position = "none",
      axis.text.y = element_markdown()
    )
  
}

# Function to plot ROC curves
plot_roc <- function(auc_data, roc_data, tag = NULL){
  
  legend_labels <- glue("{auc_data$location} (AUC = {round(auc_data$auc, 2)})")
  names(legend_labels) <- auc_data$model  # for mapping to model factor levels in plot
  
  ggplot(roc_data, aes(x = 1 - specificity, col = model)) +
    geom_step(aes(y = sensitivity), linewidth = 1.2) +
    geom_abline(linetype = "dashed") +
    scale_color_manual(
      values = lake_colors[auc_data$model],  # consistent mapping
      breaks = auc_data$model,
      labels = legend_labels
    ) +
    guides(color = guide_legend(title = NULL)) +
    labs(tag = tag) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(.99, 0.01),
      legend.justification.inside = c(1, 0),
      legend.text = element_text(size = 12),
      legend.key.spacing = unit(1, "pt"),
      legend.margin = margin(t = 0, r = 2, b = 4, l = 0, unit = "pt"),
      plot.margin = margin(t = 5, r = 5, b = 15, l = 15, unit = "pt"),
      legend.title = element_blank()
    )
}


# Read data ---------------------------------------------------------------

grab_var_imp <- read_csv("results/var-imp-grab.csv") |> 
  mutate(target = "grab")

passive_var_imp <- read_csv("results/var-imp-passive.csv") |> 
  mutate(target = "passive")

roc_grab <- read_csv("results/grab-roc-data.csv")
roc_passive <- read_csv("results/passive-roc-data.csv")  

auc_grab <- read_csv("results/grab-auc.csv") 
auc_passive <- read_csv("results/passive-auc.csv") 

full_data <- read_csv("data-cleaned/data-complete.csv") |> 
  mutate(
    location = fct_recode(location,
                          "SC" = "Shubenacadie Canal",
                          "TRR" = "Tower Road Reservoir",
                          "TCR" = "Turtle Creek Reservoir")
  )

# Organize ----------------------------------------------------------------

var_imp_data <- bind_rows(grab_var_imp, passive_var_imp) 

# Model-location lookup
location_lookup <- c(
  "full_test_kearney" = "Kearney Lake",
  "full_test_b2" = "Lake Banook 2",
  "full_test_cunard" = "Cunard Pond",
  "full_test_lc1"     = "Lake Charles 1",
  "full_test_lc2"     = "Lake Charles 2",
  "full_test_fletcher"= "Lake Fletcher",
  "full_test_sc"      = "SC",
  "full_test_inp"     = "INP",
  "full_test_trr"     = "TRR",
  "full_test_tcr"     = "TCR",
  "full_test_p" = "Penhorn Lake"
)

model_location_tbl <- tibble(
  model = names(location_lookup),
  location = unname(location_lookup)
)

# Nicer names -------------------------------------------------------------

# Vector of nicer looking names
variable_names <- c("Location" = "location",
                    "*mcyE*" = "mcye_ndaf",
                    "Cyano 16S rRNA" = "total_cyano",
                    "Temperature" = "temp_field",
                    "DO" = "do_field",
                    "Conductivity" = "cond_field",
                    "pH" = "ph_field",
                    "TDS" = "tds_field",
                    "Aluminum" = "al",
                    "Total P" = "p",
                    "Iron" = "fe",
                    "Chloride" = "chloride",
                    "NO<sub>3</sub>" = "nitrate",
                    "UV<sub>254</sub>" = "uv_254",
                    "Turbidity" = "turbidity",
                    "Colour" = "colour",
                    "DOC" = "doc",
                    "TOC" = "toc",
                    "Total N" = "tn",
                    "SO<sub>4</sub>" = "sulfate",
                    "Daily Max Temp" = "max_temp_c",
                    "Daily Min Temp" = "min_temp_c",
                    "Daily Avg Temp" =  "mean_temp_c",
                    "7d Avg Temp" =  "temp_7d_mean",
                    "3d Avg Temp" = "temp_3d_mean",
                    "7d Avg High" =  "temp_7d_avgmax",
                    "3d Avg High" =  "temp_3d_avgmax",
                    "Total Precip" = "total_precip_mm", 
                    "7d Total Precip" = "precip_7d_sum",
                    "3d Total Precip" = "precip_3d_sum",
                    "P:N" = "tn_tp_ratio",
                    "Lag Cyano 16S rRNA" = "lag_tc",
                    "Lag *mcyE*" = "lag_mcye")

# Vector of nicer target names
target_names <- c("Grab sample MC-LR" = "grab",
                  "Passive sample MC-LR" = "passive")


lake_colors <- c(
  "full_test_kearney" = "#1F77B4",
  "full_test_b2" = "#FF7F0E",
  "full_test_cunard" = "#2CA02C",
  "full_test_lc1" = "#D62728",
  "full_test_fletcher" = "#9467BD",
  "full_test_sc" = "#8C564B",
  "full_test_lc2" = "#E377C2",
  "full_test_inp" = "#7F7F7F",
  "full_test_trr" = "#BCBD22",
  "full_test_tcr" = "#17BECF",
  "full_test_p" = "#F781BF"
)

# Variable importance plots -----------------------------------------------

# Get the top n predictors for each sample type
top_n <- 15

top_grab <- var_imp_data |> 
  filter(model_type != "Tested on 2024", target == "grab") |> 
  slice_max(avg_importance, n = top_n) |> 
  mutate(
    variable = fct_reorder(variable, avg_importance),
    variable = fct_recode(variable, !!!variable_names)
  )

top_passive <- var_imp_data |> 
  filter(model_type != "Tested on 2024", target == "passive") |> 
  slice_max(avg_importance, n = top_n) |> 
  mutate(
    variable = fct_reorder(variable, avg_importance),
    variable = fct_recode(variable, !!!variable_names),
  )

# Find variables that are in the top 15 of each
shared_vars <- intersect(top_grab$variable, top_passive$variable)

# Flag them in each df
top_grab <- top_grab |> 
  mutate(
    shared = if_else(variable %in% shared_vars, "Shared predictor", "Not shared"),
    shared = factor(shared, levels = c("Shared predictor", "Not shared")))

top_passive <- top_passive |> 
  mutate(shared = if_else(variable %in% shared_vars, "Shared predictor", "Not shared"),
         shared = factor(shared, levels = c("Shared predictor", "Not shared")))


# ROC curves --------------------------------------------------------------

# Only include lakes that had at least 3 detections in the ROC plots
grab_locations <- full_data |> 
  group_by(location) |> 
  summarise(n = sum(grab_mclr > grab_mclr_lod, na.rm = TRUE)) |> 
  filter(n > 2) |> 
  inner_join(model_location_tbl, by = "location")

passive_locations <- full_data |> 
  group_by(location) |> 
  summarise(n = sum(passive_mclr > passive_mclr_lod, na.rm = TRUE)) |> 
  filter(n > 2) |> 
  inner_join(model_location_tbl, by = "location")

# Combine plots -----------------------------------------------------------

grab_importance_plot <- 
  var_imp_plotting(top_grab) + 
  labs(alpha = NULL, tag = "b)") + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.margin =  margin(b = -10, unit="pt"),
  )

passive_importance_plot <- var_imp_plotting(top_passive) + labs(tag = "d)")

grab_roc_plot <- plot_roc(
  auc_data = auc_grab |> filter(model != "full_test_2024", location %in% grab_locations$location), 
  roc_data = roc_grab |> filter(model %in% grab_locations$model), tag = "a)") + 
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
    )

passive_roc_plot <- plot_roc(
  auc_data = auc_passive |> filter(location %in% passive_locations$location),
  roc_data = roc_passive |> filter(model %in% passive_locations$model),
  tag = "c)") +
  guides(color = guide_legend(ncol = 2, title = NULL)) 


((grab_roc_plot / passive_roc_plot) | (grab_importance_plot / passive_importance_plot)) +  plot_layout(ncol = 2, widths = c(6, 2.5))  &
  theme(plot.tag = element_text(face = "bold"))


ggsave("figures/roc-var-importance-passive-grab.png", width = 15, height = 12, dpi = 300)


# Plot the models trained and tested on year sets -------------------------

# Get the top n predictors for each sample type

top_grab_year <- var_imp_data |> 
  filter(model_type == "Tested on 2024", target == "grab") |> 
  slice_max(avg_importance, n = top_n) |> 
  mutate(
    variable = fct_reorder(variable, avg_importance),
    variable = fct_recode(variable, !!!variable_names)
  )

top_passive_year <- var_imp_data |> 
  filter(model_type == "Tested on 2024", target == "passive") |> 
  slice_max(avg_importance, n = top_n) |> 
  mutate(
    variable = fct_reorder(variable, avg_importance),
    variable = fct_recode(variable, !!!variable_names),
  )

# Find variables that are in the top 15 of each
shared_vars_yr <- intersect(top_grab_year$variable, top_passive_year$variable)

# Flag them in each df
top_grab_year <- top_grab_year |> 
  mutate(
    shared = if_else(variable %in% shared_vars_yr, "Shared predictor", "Not shared"),
    shared = factor(shared, levels = c("Shared predictor", "Not shared")))

top_passive_year <- top_passive_year |> 
  mutate(shared = if_else(variable %in% shared_vars_yr, "Shared predictor", "Not shared"),
         shared = factor(shared, levels = c("Shared predictor", "Not shared")))


# ROC curves for year models ----------------------------------------------

grab_importance_plot_yr <- 
  var_imp_plotting(top_grab_year) + 
  labs(alpha = NULL, tag = "b)") + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.margin =  margin(b = -10, unit="pt"),
  )

passive_importance_plot_yr <- var_imp_plotting(top_passive_year) + labs(tag = "d)")

roc_year_grab <- roc_grab |>
  filter(model_type == "Tested on 2024") |> 
  ggplot(aes(x = 1-specificity)) +
  pammtools::geom_stepribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = "grey97") +
  geom_step(aes(y = sensitivity, col = model), linewidth = 1.2) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(
    # breaks = c("full", "subset"),
    values = "red",
    labels = glue::glue("2024 (AUC = {auc_grab[auc_grab$model == 'full_test_2024', ]$auc})")) +
  guides(
    fill = "none",
    color = guide_legend(
      title = NULL,
      position = "inside")) +
  labs(tag = "a)") +
  theme(
    legend.position.inside = c(.99, 0.01),
    legend.justification.inside = c(1,0))

roc_year_passive <- roc_passive |>
  filter(model_type == "Tested on 2024") |> 
  ggplot(aes(x = 1-specificity)) +
  pammtools::geom_stepribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = "grey97") +
  geom_step(aes(y = sensitivity, col = model), linewidth = 1.2) +
  geom_abline(linetype = "dashed") +
  scale_color_manual(
    # breaks = c("full", "subset"),
    values = "red",
    labels = glue::glue("2024 (AUC = {auc_passive[auc_passive$model == 'full_test_2024', ]$auc})")) +
  guides(
    fill = "none",
    color = guide_legend(
      title = NULL,
      position = "inside")) +
  labs(tag = "c)") +
  theme(
    legend.position.inside = c(.99, 0.01),
    legend.justification.inside = c(1,0))

((roc_year_grab / roc_year_passive) | (grab_importance_plot_yr / passive_importance_plot_yr)) +  plot_layout(ncol = 2, widths = c(6, 2.5))  &
  theme(plot.tag = element_text(face = "bold"))


ggsave("figures/si-figures/roc-var-importance-passive-grab-year.png", width = 15, height = 12, dpi = 300)
