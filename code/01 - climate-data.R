# Looking at climate data to see if there are any weather trends that might help explain 
# MC-LR detections


# Setup -------------------------------------------------------------------

library("tidyverse")
library("janitor")

# Read NS Environment files -----------------------------------------------

climate_files <- list.files("data-raw/climate-data", full.names = TRUE) %>%
  tibble() |> 
  filter(str_detect(., "en_climate"))

climate_data <- map(climate_files, read_csv) |>
  list_rbind() |> 
  clean_names() |> 
  mutate(month = as.numeric(month)) |> 
  filter(month %in% 4:12) |> 
  select(
    -c(longitude_x, latitude_y, climate_id,
       year, month, day, data_quality),
    -contains("flag"), -contains("gust"), -c(snow_on_grnd_cm, total_rain_mm, total_snow_cm)) |> 
  pivot_longer(
    -c(station_name, date_time),
    names_to = "param",
    values_to = "value"
  ) |> 
  group_by(station_name, date_time, param) |> 
  summarise(value = mean(value, na.rm = TRUE)) |> 
  ungroup() |> 
  pivot_wider(
    id_cols = c(station_name, date_time),
    names_from = param,
    values_from = value
  ) |> 
  mutate(region = if_else(str_detect(station_name, "MONCTON"), "nb", "ns")) 

write_csv(climate_data, "data-raw/climate-data-compiled.csv")
