# plots that explore the relationship between nutrients and MC-LR detections

library(tidyverse)
library(assertr)
library(ggtext)

theme_set(theme_bw(base_size = 16))

detection_color <- "firebrick"
grab_col <- "#f7d13d"
passive_col <- "#4a0c6b"

# Named vectors for plotting ----------------------------------------------

# Params to plot
these_params <- c("p_ugl", "tn_ugl", "colour", "toc_mgl", "tn_tp_ratio")

# Tidy labels
facet_labs <- list(
  analyte = c(
    # chloride_mgl = "Chloride (mgL<sup>-1</sup>)",
    p_ugl = "Total Phosphorus (&mu;gL<sup>-1</sup>)",
    tn_ugl = "Total Nitrogen (&mu;gL<sup>-1</sup>)",
    colour = "True Colour (ptCo)",
    # phosphate_ugl = "Phosphate (&mu;gL<sup>-1</sup> as P)",
    toc_mgl = "Total Organic Carbon (mgL<sup>-1</sup>)",
    tn_tp_ratio = "TN:TP"
    # doc_mgl = "Dissolved Organic Carbon (mgL<sup>-1</sup>)"
    )
  )

# Swap names and values for fct_recode()
param_labels <- setNames(names(facet_labs$analyte), facet_labs$analyte)
# Add a line break 
names(param_labels) <- names(param_labels) |> 
  str_replace( "(\\(.*\\))", "<br>\\1")

# Read data ---------------------------------------------------------------

complete_data <- read_csv(here::here("data-cleaned/data-complete.csv")) |> 
  mutate(
    year = as.factor(year),
    chloride_mgl = chloride / 1000,
    tn_ugl = tn * 1000,
    p_ugl = p,
    toc_mgl = toc,
    doc_mgl = doc,
    total_cyano = log10(total_cyano + 1e-2),
    mcye_ndaf = log10(mcye_ndaf + 1e-2),
    # create nutrient ratio parameters (TN and TOC are in mg/L, TP is µg/L)
    tn_tp_ratio = (tn / 14) / ((p / 1000) / 31),
    tc_tn_ratio = (toc / 12) /  (tn / 14),
    tc_tp_ratio = (toc / 12) / ((p / 1000) / 31)
    )

# Locations with most MC-LR detections and max grab concentrations --------

most_detections_df <- complete_data |> 
  group_by(location, sample_date) |>
  mutate(n = n(), .before= everything()) |> 
  verify(n == 1) |>
  ungroup() |> 
  select(-n) |> 
  # remove the missing samples
  filter(!is.na(passive_mclr)) |> 
  group_by(location) |> 
  summarise(
    n = n(),
    n_detections = sum(passive_mclr > passive_mclr_lod, na.rm = TRUE),
    det_rate = paste(n_detections, "/", n, sep = " ")) |> 
  arrange(desc(n_detections))  |> 
  mutate(
    location = fct_recode(location,
                        "SC" = "Shubenacadie Canal",
                        "TRR" = "Tower Road Reservoir",
                        "TCR" = "Turtle Creek Reservoir"),
  label = paste(location,"<br>" ,"(", det_rate,")", sep = ""))
  
most_detections <- most_detections_df |> pull(label)

# Correlation between detection rate and median nutrient concentrations----

spearman_correlations <- complete_data |>
  filter(!is.na(passive_mclr)) |> 
  select(location, all_of(these_params), contains("passive_mclr")) |> 
  # indicator column for MC-LR detections
  mutate(detected = passive_mclr > passive_mclr_lod) |> 
  select(-contains("mclr")) |> 
  group_by(location) |> 
  summarise(
    across(p_ugl:tn_tp_ratio, ~ median(.x, na.rm = TRUE)),
    n = n(),
    detection_rate = mean(detected)
  ) |> 
  pivot_longer(-c(location, n, detection_rate)) |> 
  group_by(name) |> 
  nest() |> 
  ungroup() |> 
  mutate(
    spearman = map(data, ~ cor(x = .x$value, y = .x$detection_rate, method = "spearman"))
  ) |> 
  select(name, spearman) |> 
  unnest(spearman) |> 
  mutate(
    spearman = round(spearman, digits = 2),
    name = factor(name, levels = 
                    c("chloride_mgl", "p_ugl", "tn_ugl", "colour", "toc_mgl", "doc_mgl", "tn_tp_ratio")))

# Faceted jitter plot -----------------------------------------------------

# generate plot
plot_data <- complete_data |>
  select(location, year, all_of(these_params), contains("passive_mclr")) |> 
  # indicator column for MC-LR detections
  mutate(
    detected = if_else(passive_mclr > passive_mclr_lod, "Detected", "Not detected"),
    detected = factor(detected, levels = c("Not detected", "Detected", NA)),
    location = fct_recode(location,
                          "SC" = "Shubenacadie Canal",
                          "TRR" = "Tower Road Reservoir",
                          "TCR" = "Turtle Creek Reservoir")) |> 
  pivot_longer(-c(location, year, detected)) |> 
  filter(
    name %in% these_params,
    !is.na(detected),
    case_when(
      name == "p_ugl" & value > 250 ~ FALSE,
      name == "tn_ugl" & value > 1000 ~ FALSE,
      # name == "phosphate_ugl" & value > 50 ~ FALSE,
      TRUE ~ TRUE
    )
  ) |> 
  left_join(most_detections_df, by = "location") |>
  # left_join(max_concentrations_df, by = "location") |> 
  mutate(
    # order the locations by the ordering found above
    label = factor(label, levels = rev(most_detections)),
    # label = factor(label, levels = rev(max_concentrations)),
    name = factor(name, levels = 
                    c("chloride_mgl", "p_ugl", "tn_ugl", "colour", "toc_mgl", "doc_mgl", "tn_tp_ratio")),
    n_detections = factor(n_detections)
    # max = factor(max)
  ) 
  

nutrients_toxins_plot <- ggplot(plot_data, aes(label, value)) +
  facet_wrap(vars(name), nrow = 1, scales = "free_x",
             labeller = labeller(name = facet_labs$analyte)) +
  geom_jitter(aes(color = detected,  alpha = detected, size = detected), 
              show.legend = c(color = TRUE, alpha = FALSE, size = FALSE),
              width = 0.1) +
  # geom_boxplot(alpha = 0.1, width = 0.5) +
  coord_flip() +
  scale_color_manual(values = c("grey", detection_color)) +
  scale_alpha_manual(values = c(0.4, 0.9)) +
  scale_size_manual(values = c(2,3)) +
  guides(
    color = guide_legend(
      title = "MC-LR",
      position = "bottom",
      override.aes = list(size = 3)
      )
    ) +
  ggh4x::facetted_pos_scales(
    y = list(
      name == "toc_mgl" ~ scale_y_continuous(limits = c(0,15)),
      name == "p_ugl" ~ scale_y_continuous(expand = c(0.1, 0.2))
    )
  ) +
  geom_richtext(
    data = spearman_correlations,
    aes(x = Inf, y = Inf,
        label = paste("*&rho;* = ",spearman, sep = ""),
        hjust = 1, vjust =1,
        fontface = "bold")
  ) +
  theme(
    axis.text.y = element_markdown(hjust = 1),
    strip.text.x = element_markdown()
  ) +
  labs(
    x = NULL,
    y = NULL
  )

nutrients_toxins_plot

ggsave("figures/nutrients-and-toxins.png", width = 14, height = 8)


