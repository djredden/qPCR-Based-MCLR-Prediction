# Figure that shows the absolute range of MC-LR in grab samples and the 
# passive sample MC-LR detections while giving a clear representation of the 
# number of passive samples at each location (important to note the lost samplers)

# Setup -------------------------------------------------------------------

library("tidyverse")
library("ggtext")
library("patchwork")

theme_set(theme_bw(base_size = 16))

detection_color <- "firebrick"
grab_col <- "#51c016"
passive_col <- "#4a0c6b"
year_colors <- c("#E69F00","#56B4E9", "#009E73")

# Read --------------------------------------------------------------------

complete_data <- read_csv("data-cleaned/data-complete.csv") |> 
  mutate(year = as.factor(year)) 

# Do passive samples precede grab samples? --------------------------------

# Get toxin data in long format
passive_grab <- complete_data |>
  select(location, sample_date, contains("grab"), contains("passive")) |>
  rename(
    grab_value = grab_mclr,
    passive_value = passive_mclr,
    grab_lod = grab_mclr_lod,
    passive_lod = passive_mclr_lod) |>
  pivot_longer(-c(location, sample_date),
               names_to = c("method", ".value"),
               names_sep = "_")

# Filter grab detections
grab_detections <- passive_grab  |> 
  filter(method == "grab", value > lod) |> 
  select(location, sample_date) |> 
  mutate(grab_detection = TRUE)

# Check if there were passive detections the same week or in the previous two weeks
paired_detections <- passive_grab |> 
  filter(method == "passive") |> 
  left_join(grab_detections, by = c("location", "sample_date")) |> 
  group_by(location, year = year(sample_date)) |> 
  mutate(
    passive_detection = value > lod,
    # passive_detection = replace_na(passive_detection, FALSE),
    one_week = lag(passive_detection),
    two_weeks = lag(passive_detection, 2)
  ) |> 
  ungroup() |> 
  filter(grab_detection) |> 
  group_by(location, year) |> 
  summarise(
    n_detects = sum(grab_detection, na.rm = TRUE),
    same_week = sum(passive_detection, na.rm = TRUE),
    prev_week = sum(one_week, na.rm = TRUE),
    two_weeks = sum(two_weeks,  na.rm = TRUE),
    .groups = "drop") 

# Plot grab and passive time series together ------------------------------

# Specify date breaks
axis_dates <- as.Date(c("2022-05-01", "2022-08-01", "2022-11-01",
                          "2023-05-01", "2023-08-01", "2023-11-01",
                          "2024-05-01", "2024-08-01", "2024-11-01"
                          ))

year_labels <- data.frame(
  sample_date = as.Date(c("2022-08-01", "2023-08-01", "2024-08-01")),
  value = 0,  # or use slightly negative if your y-axis starts above zero
  label = c("2022", "2023", "2024")
)
  
passive_grab_plot <- passive_grab |>
  # remove missing data and censored grab data
  filter(
    !is.na(value), 
    case_when(
      method == "grab" & value <= lod ~ FALSE,
      TRUE ~ TRUE)
    ) |> 
  mutate(
    cens = if_else(value <= lod, "&lt; LOD", "&ge; LOD"),
    cens = factor(cens, levels = c("&ge; LOD", "&lt; LOD")),
    method = fct_recode(method,
                        "Grab MC-LR (&mu;g L<sup>-1</sup>)" = "grab",
                        "Passive MC-LR (&mu;g g<sup>-1</sup>)" = "passive")) |> 
  ggplot(aes(sample_date, value, alpha = cens, col = method)) + 
  facet_wrap(vars(location), ncol = 2) +
  geom_vline(
    xintercept = as.numeric(axis_dates),  # convert to numeric to match scale_x_date
    color = "grey90"
  ) +
  geom_point(size = 3) +
  geom_vline(
    data = . %>% filter(method == "Grab MC-LR (&mu;g L<sup>-1</sup>)"),
    aes(xintercept = sample_date),
    linetype = "42", color = "grey50"
  ) +
  scale_alpha_manual(values = c(0.9, 0.15)) +
  scale_x_date(
    breaks = axis_dates,
    labels = c("May",  "Aug <br> 2022",  "Nov",
               "May",  "Aug <br> 2023",  "Nov",
               "May",  "Aug <br> 2024",  "Nov")
  ) +
  scale_color_manual(values = c(grab_col, passive_col)) +
  guides(
    alpha = guide_legend(
      title = element_blank(),
      position = "bottom", 
      order = 2,
      theme(
        legend.title = element_text(size = 14),
        legend.margin =  margin(t = -15, b= 0, unit="pt"),
        legend.justification.bottom = "center")
    ),
    color = guide_legend(
      title = element_blank(), 
      position = "bottom",
      order = 1,
      override.aes = list(size = 3),
      theme(
        legend.text = element_markdown(margin = margin(l = 0)),
        legend.margin =  margin(t= -10, unit="pt"),
        legend.justification.bottom = "center")
    )) +
  labs(
    x = NULL,
    y = "MC-LR Concentration") +
  theme(
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(0, "lines"),
    legend.box = "vertical",
    strip.text.y = element_text(angle = 0),
    strip.clip = "off",
    panel.border = element_rect(linewidth = 0.5),
    strip.background = element_rect(linewidth = 0.5, fill = "lightgrey"),
    legend.text = element_markdown(),
    axis.text.x = element_markdown()
    # axis.text.x = element_markdown(angle = 45, vjust = 0.75, hjust = 0.75)
  ) 

passive_grab_plot

# Plot grab sample concentrations -----------------------------------------

# Some additional wrangling
these_locations <- complete_data |>
  filter(grab_mclr > grab_mclr_lod) |> 
  distinct(location) |> # locations where MC-LR was detected in grabs
  pull(location)

# Plot the grab concentrations in each year -------------------------------

grab_plot <- complete_data |> 
  mutate(
    year = year(sample_date),
    year = as.factor(year),
    grab_cens = if_else(is.na(grab_mclr) | grab_mclr > grab_mclr_lod, "none", "left"),
    month_date = as_date(paste("2022", month(sample_date), day(sample_date), sep = "-"))) |> 
  filter(
    location %in% these_locations,
    grab_cens == "none") |> 
  ggplot(aes(month_date, grab_mclr, shape = location, color = year)) +
  geom_point(size = 3.5, stroke = 2, alpha = 0.6) +
  scale_shape_manual(values = c(0,1,2,3,4,15,16,17,8,9)) +
  scale_x_date(
    date_breaks = "month",
    date_labels = "%b") +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = year_colors) +
  guides(
    col = guide_legend(
      title = NULL, 
      position = "inside", 
      theme(
        legend.background = element_rect(color = "black")
        )
      ),
    shape = guide_legend(
      position = "bottom",
      ncol = 3,
      override.aes = list(size = 3, stroke = 1)
      )
    ) +
  labs(
    y = "Grab Sample <br> MC-LR (&mu;g L<sup>-1</sup>)",
    x = NULL,
    shape = NULL,
    col = NULL
  ) +
  theme(
    axis.title.y = element_markdown(),
    legend.position.inside = c(0.1, 0.9),
    legend.key.spacing.y = unit(2, "pt"))
  
ggsave("figures/si-figures/grab-mclr.png", grab_plot,
       dev = "png", dpi = 600,
       width = 8, height = 8)


# free(grab_plot, side = "t") +  plot_spacer() + (passive_grab_plot) + plot_layout(ncol = 3, widths = c(2,0.1, 6)) &
#   theme(plot.tag = element_text(face = "bold"))


ggsave("figures/passive-mclr.png", passive_grab_plot,
       dev = "png", dpi = 600,
       width = 16, height = 8)


