
# Impute missing data from the lake monitoring program using random forest.

library(tidyverse)
library(ggtext)
library(missRanger)

# Read data ---------------------------------------------------------------

data_raw <- read_csv("data-cleaned/data-complete.csv") # data contains only predictors with <40% obs missing

# Select those variables to be used for analysis ---------------------------

toxin_vars <- c("grab_mclr", "passive_mclr", "spe_mclr") # toxin variables

model_params <- c("location", "sample_date", "temp_field", "do_field", "cond_field", "ph_field", 
                  "tds_field" , "al", "p", "fe", "chloride", "nitrate", "sulfate", "uv_254", 
                  "turbidity", "colour", "doc", "toc", "tn", "cyra", "mcye_ndaf", "sxta", 
                  "total_cyano", "max_temp_c", "min_temp_c", "mean_temp_c", "temp_7d_mean",
                  "temp_3d_mean", "temp_7d_avgmax", "temp_3d_avgmax", "total_precip_mm", 
                  "precip_7d_sum", "precip_3d_sum" )

censoring_params <- c("ul_al", "rdl_al", "ul_p", "rdl_p", "ul_fe", "rdl_fe", "lod_cyra", "lod_mcye_ndaf", 
                      "lod_sxta", "lod_total_cyano", "grab_mclr_lod", "passive_mclr_lod", "spe_mclr_lod", 
                       "rdl_chloride",  "rdl_nitrate", "rdl_sulfate")

# Subset the data to those params
data_to_impute <- data_raw |> 
  mutate(location = factor(location)) |> 
  select(all_of(model_params)) # don't want to impute the responses

# MissRanger imputation ---------------------------------------------------

# Set up the imputation
rf_imputation <- vector(mode = "list") # initialize list
n_datasets <- 25 # number of imputed datasets to create

# impute and join the MCLR data back to each imputed df
for(i in 1:n_datasets) {
  rf_imputation[[i]] <- missRanger(data_to_impute, num.trees = 100, pmm.k = 5)
  # add the toxin data and censoring limits to each dataset
  rf_imputation[[i]] <- bind_cols(
    data_raw |> select(all_of(toxin_vars), all_of(censoring_params)),
    rf_imputation[[i]])
}

saveRDS(rf_imputation, "model-fits/miss_ranger_imputation.rds")
rf_imputation <- readRDS("model-fits/miss_ranger_imputation.rds")

# Plot data distributions before and after imputation ---------------------

nice_names <- tibble(
  original_name = names(data_to_impute),
  nice_name = c("Location",
                "Sample Date",
                "Water Temperature (&deg;C)", 
                "Dissolved Oxygen (mgL<sup>-1</sup>)", 
                "Conductivity (&mu;S cm<sup>-1</sup>)",
                "pH", 
                "TDS (mg L<sup>-1</sup>)", 
                "Aluminum (&mu;gL<sup>-1</sup>)",
                "Total Phosphorus (&mu;gL<sup>-1</sup>)",
                "Iron (&mu;gL<sup>-1</sup>)",
                "Chloride (mgL<sup>-1</sup>)",
                "Nitrate (&mu;gL<sup>-1</sup> as N)",
                "Sulfate (mgL<sup>-1</sup>)",
                "UV<sub>254</sub> (cm<sup>-1</sup>)",
                "Turbidity (NTU)",
                "True colour (ptCo)",
                "Dissolved Organic Carbon (mg L<sup>-1</sup>)",
                "Total Organic Carbon (mg L<sup>-1</sup>)",
                "Total Nitrogen (&mu;g L<sup>-1</sup>)",
                "*cyr*A (GC mL<sup>-1</sup>)",
                "*mcy*E (GC mL<sup>-1</sup>)",
                "*sxt*A (GC mL<sup>-1</sup>)",
                "Total Cyanobacteria (GC mL<sup>-1</sup>)",
                "Max air temp (&deg;C)",
                "Min air temp (&deg;C)",
                "Avg air temp (&deg;C)",
                "7d avg air temp (&deg;C)",
                "3d avg air temp (&deg;C)",
                "7d avg max air temp (&deg;C)",
                "3d avg max air temp (&deg;C)",
                "Total precipitation (mm)",
                "7d total precipitation (mm)",
                "3d total precipitation (mm)"))

# Parameter names to plot
plot_vars <- data_to_impute |> 
  names() %>% 
  tibble("params" = .) |> 
  filter(
    !params %in% c("location", "sample_date")
  ) |> 
  pull()


imputed_data <- list_rbind(rf_imputation)

distribution_plot_function <- function(param) {
  
  x_axis_title <- nice_names |> 
    filter(original_name == param) |> 
    pull(nice_name)
  
  ggplot() +
    facet_wrap(vars(location), scales = "free") +
    
    # Original data
    geom_density(data = data_to_impute, 
                 aes(x = .data[[param]], fill = "Original", group = location), 
                 alpha = 0.6) +
    
    # Imputed data
    geom_density(data = imputed_data, 
                 aes(x = .data[[param]], fill = "Imputed", group = location), 
                 alpha = 0.6) +
    
    labs(x = x_axis_title, 
         y = NULL,
         fill = NULL) +
    theme_bw() +
    theme(axis.title.x = element_markdown())
}

# create plots
distibution_plots <- map(plot_vars, distribution_plot_function)

names(distibution_plots) <- plot_vars

distibution_plots %>%
  names(.) |> 
  walk(~ ggsave(paste0("extra-figures/imputation-distributions/", ., ".png"), 
                distibution_plots[[.]],  
                height = 6, width = 14))

# Plot mcyE and total cyano separately to allow log10 transformation

mcye_dist <- ggplot() +
  facet_wrap(vars(location), scales = "free") +
  
  # Original data
  geom_density(data = data_to_impute, 
               aes(x = log10(.data[["mcye_ndaf"]]), fill = "Original", group = location), 
               alpha = 0.6) +
  
  # Imputed data
  geom_density(data = imputed_data, 
               aes(x = log10(.data[["mcye_ndaf"]]), fill = "Imputed", group = location), 
               alpha = 0.6) +
  
  labs(x = "log<sub>10</sub> *mcyE* (GC mL<sup>-1</sup>)", 
       y = NULL,
       fill = NULL) +
  theme_bw() +
  theme(axis.title.x = element_markdown())

ggsave("extra-figures/imputation-distributions/mcye_ndaf.png", 
       mcye_dist, 
       height = 6, width = 14)


tc_dist <- ggplot() +
  facet_wrap(vars(location), scales = "free") +
  
  # Original data
  geom_density(data = data_to_impute, 
               aes(x = log10(.data[["total_cyano"]]), fill = "Original", group = location), 
               alpha = 0.6) +
  
  # Imputed data
  geom_density(data = imputed_data, 
               aes(x = log10(.data[["total_cyano"]]), fill = "Imputed", group = location), 
               alpha = 0.6) +
  
  labs(x = "log<sub>10</sub> Total Cyanobacteria (GC mL<sup>-1</sup>)", 
       y = NULL,
       fill = NULL) +
  theme_bw() +
  theme(axis.title.x = element_markdown())

ggsave("extra-figures/imputation-distributions/total_cyano.png", 
       tc_dist, 
       height = 6, width = 14)
