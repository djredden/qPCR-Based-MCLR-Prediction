# Notes: 

# 1. MC-LR values — grab, passive, and SPE — in qpcr-algal-toxin-data.csv have been replaced by (censored at)
#    their respective LODs where they were <LOD
# 2. In some cases solid phase extraction (SPE) replaced typical grab samples. In those case, the grab MC-LR value was replaced
#    by the SPE value but censored at the grab LOD
# 3. Missing MC-LR values represent true missing samples
# 4. Max values of toxins have been used wherever there were multiple samples collected in a week
# 5. Max values of qPCR results have been used wherever there were multiple samples in a week

library("tidyverse")
source("code/functions.R")

# Read data ---------------------------------------------------------------

data_raw <- read_csv("data-raw/qpcr-algal-toxin-data.csv") 
spe_raw <- read_csv("data-raw/spe-data.csv")
climate_data_raw <- read_csv("data-raw/climate-data-compiled.csv")

# Create nicer names ------------------------------------------------------

data <- data_raw |> 
  mutate(region = if_else(location %in% c("INP", "Tower Road", "Turtle Creek"), "nb", "ns"),
         location = fct_recode(location,
                               "Lake Banook 1" = "Banook 1",
                               "Lake Banook 2" = "Banook 2",
                               "Shubenacadie Canal" = "Fairbanks",
                               "Penhorn Lake" = "Penhorn",
                               "Lake Charles 1" = "Shubie 1",
                               "Lake Charles 2" = "Shubie 2",
                               "Lake Fletcher" = "Collins Park",
                               "Kearney Lake" = "Kearney",
                               "Cunard Pond" = "Cunard",
                               "Tower Road Reservoir" = "Tower Road",
                               "Turtle Creek Reservoir" = "Turtle Creek",
                               "INP" = "INP"))

spe <- spe_raw |>
  mutate(location = fct_recode(location,
                               "Lake Charles 2" = "Shubie 2",
                               "Lake Fletcher" = "Collins Park",
                               "Cunard Pond" = "Cunard",
                               "INP" = "INP"))

# Replace missing grab samples with SPE and censor at grab LOD ---------

complete_data <- data |> 
  mutate(grab_imputed = if_else(is.na(grab_mclr), spe_grab_mclr, grab_mclr)) |> 
  select(-grab_mclr) |> 
  rename(
    spe_mclr = spe_grab_mclr,
    spe_mclr_lod = spe_grab_mclr_lod,
    grab_mclr = grab_imputed
  ) |> 
  relocate(c(grab_mclr, passive_mclr, spe_mclr), .after = sample_time) 


# Assess missingness ------------------------------------------------------

# Parameters missing more than 40% of the time are excluded from analyses

# Proportion of missing values for each parameter
prop_missing <- complete_data |>
  select(-c(id, location, sample_date, sample_time, region),
         -contains("rdl"),
         -contains("ul_"), 
         -contains("lod"),
         -sample) |>
  pivot_longer(everything()) |> 
  group_by(name) |> 
  summarise(prop_missing = mean(is.na(value))) |> 
  arrange(desc(prop_missing))

# Parameters missing more than 40% of observations
hi_missing <- prop_missing |> 
  filter(
    prop_missing > 0.4,
    !str_detect(name, "spe")) |> # ensure SPE data is retained
  pull(name) 

# Organize climate data ---------------------------------------------------

climate_data <- climate_data_raw |> 
  mutate(
    precip_7d_sum = zoo::rollapply(total_precip_mm, 7, sum, align = "right", fill = NA, na.rm = TRUE),
    precip_3d_sum = zoo::rollapply(total_precip_mm, 3, sum, align = "right", fill = NA, na.rm = TRUE),
    temp_7d_mean = zoo::rollapply(mean_temp_c, 7, mean, align = "right", fill = NA, na.rm = TRUE),
    temp_3d_mean = zoo::rollapply(mean_temp_c, 3, mean, align = "right", fill = NA, na.rm = TRUE),
    temp_7d_avgmax = zoo::rollapply(max_temp_c, 7, mean, align = "right", fill = NA, na.rm = TRUE),
    temp_3d_avgmax = zoo::rollapply(max_temp_c, 3, mean, align = "right", fill = NA, na.rm = TRUE)
  ) |> 
  select(date_time, region, contains("temp"), contains("precip")) 
  

# Combine climate and wq data ---------------------------------------------

# Final dataset -----------------------------------------------------------

final_data <- complete_data |> 
  left_join(climate_data, by = c("sample_date" = "date_time", "region")) |> 
  # remove the params with high missingness and dissolved nitrogen (not typically used) and their det limits
  select(!any_of(contains(hi_missing)), -dn, -region) |> # 
  relocate(contains("ul_"), .after = everything()) |> 
  relocate(contains("lod"), .after = everything()) |> 
  relocate(contains("rdl"), .after = everything())

visdat::vis_miss(final_data)

# Write -------------------------------------------------------------------

# locations for this project:
study_locations <- c("Lake Banook 2", "Lake Fletcher",  "Cunard Pond", "Shubenacadie Canal", "INP" ,
                     "Kearney Lake", "Lake Charles 1", "Lake Charles 2", "Tower Road Reservoir",
                     "Turtle Creek Reservoir", "Penhorn Lake", "Lake Banook 1")
final_data |> 
  filter(location %in% study_locations) %>%
  write_csv(., "data-cleaned/data-complete.csv")

