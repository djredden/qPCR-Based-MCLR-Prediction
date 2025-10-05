# setup -------------------------------------------------------------------

library("tidyverse")
library("ggtext")
library("assertr")
library("patchwork")

theme_set(theme_bw(base_size = 16))

detection_color <- "firebrick"
grab_col <- "#51c016"
passive_col <- "#4a0c6b"
year_colors <- c("#E69F00","#56B4E9", "#009E73")

# Read --------------------------------------------------------------------

complete_data <- read_csv("data-cleaned/data-complete.csv") |> 
  mutate(year = as.factor(year)) 

# Plot qPCR results -------------------------------------------------------

axis_dates <- as.Date(c("2022-05-01", "2022-08-01", "2022-11-01",
                        "2023-05-01", "2023-08-01", "2023-11-01",
                        "2024-05-01", "2024-08-01", "2024-11-01"
))

year_labels <- data.frame(
  sample_date = as.Date(c("2022-08-01", "2023-08-01", "2024-08-01")),
  value = 0,  # or use slightly negative if your y-axis starts above zero
  label = c("2022", "2023", "2024")
)


complete_data |> 
  select(location, sample_date, contains("total_cyano"), contains("mcye_ndaf")) |> 
  rename(
    tc_value = total_cyano,
    tc_lod = lod_total_cyano,
    mcye_value = mcye_ndaf,
    mcye_lod = lod_mcye_ndaf
  ) |> 
  pivot_longer(
    -c(location, sample_date),
    names_to = c("target", ".value"),
    names_sep = "_"
    ) |> 
  filter(
    !is.na(value),
    value < 1e7) |> #remove one outlier from Penhorn
  mutate(
    detected = if_else(value > lod, "&ge; LOD", "< LOD"),
    detected = factor(detected, levels = c("&ge; LOD", "< LOD")),
    target = if_else(target == "tc", "Total cyanobacteria", "*mcyE*"),
    target = factor(target, levels = c("Total cyanobacteria", "*mcyE*"))
    ) |> 
  ggplot(aes(sample_date, log10(value + 1e-1), col = target, alpha = detected)) +
  facet_wrap(vars(location)) +
  geom_point(size = 2.5) +
  scale_alpha_manual(
    values = c(1, 0.2)) +
  ggsci::scale_color_jco() +
  scale_y_continuous(
    breaks = c(0,2,4,6),
    labels = scales::math_format()) +
  scale_x_date(
    breaks = axis_dates,
    labels = c("May",  "Aug <br> 2022",  "Nov",
               "May",  "Aug <br> 2023",  "Nov",
               "May",  "Aug <br> 2024",  "Nov")
  ) +
  labs(
    x = NULL,
    y = "Gene concentration (GC mL<sup>-1</sup>)",
    col = NULL,
    alpha = NULL
  ) +
  guides(
    color = guide_legend(order = 1),
    alpha = guide_legend(order = 2)
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.text = element_markdown(),
    legend.margin =  margin(b= -10, unit="pt"),
    plot.margin = margin(l = 5, r = 15),
    )

ggsave("figures/qpcr.png",
       dev = "png", dpi = 300,
       width = 18, height = 10)


# Kruskal Wallis test -----------------------------------------------------

# kruskall wallis
kw_results <- complete_data |> 
  mutate(location = factor(location)) |> 
  select(location, sample_date, mcye_ndaf, total_cyano) |> 
  pivot_longer(-c(location, sample_date)) |> 
  group_by(name) |> 
  nest() |> 
  ungroup() |> 
  mutate(
    kw = map(data, ~ kruskal.test(value ~ location, data = .)),
    results = map(kw, broom::glance)
  ) |> 
  select(name, results) |> 
  unnest(results) |> 
  mutate(
    p.value = round(p.value, 2),
    p_value_tidy = if_else(p.value < 0.05, "<0.05", as.character(p.value))
  )

# Plot boxplots with jitter and KW results --------------------------------

qpcr_boxplots <- complete_data |> 
  select(location, sample_date, contains("total_cyano"), contains("mcye_ndaf")) |> 
  rename(
    tc_value = total_cyano,
    tc_lod = lod_total_cyano,
    mcye_value = mcye_ndaf,
    mcye_lod = lod_mcye_ndaf
  ) |> 
  pivot_longer(
    -c(location, sample_date),
    names_to = c("target", ".value"),
    names_sep = "_"
  ) |> 
  filter(
    !is.na(value),
    value < 1e7) |> #remove one outlier from Penhorn
  mutate(
    detected = if_else(value > lod, "&ge; LOD", "< LOD"),
    detected = factor(detected, levels = c("&ge; LOD", "< LOD")),
    target = if_else(target == "tc", "Total cyanobacteria", "*mcyE*"),
    target = factor(target, levels = c("Total cyanobacteria", "*mcyE*")),
    year = factor(year(sample_date))
  ) |> 
  ggplot(aes(log10(value + 1e-1), location)) +
  facet_wrap(vars(target)) +
  geom_jitter(aes(color = year), height = 0.1, alpha = 0.6, size = 3) +
  geom_boxplot(alpha = 0) +
  expand_limits(x = 0) +
  scale_x_continuous(
    breaks = c(0,2,4,6),
    labels = scales::math_format()) +
  scale_color_manual(values = year_colors) +
  geom_richtext(
    data = kw_results,
    aes(x = Inf, y = Inf,
        label = paste("Kruskal-Wallis *p* ", p_value_tidy, sep = ""),
        hjust = 1, vjust =1)
  ) +
  guides(
    color = guide_legend(
      title = element_blank(),
      position = "top",
      override.aes = list(alpha = 1),
      theme(
        legend.margin =  margin(b = -10, unit="pt"),
        legend.justification.bottom = "center")
    )
  ) +
  labs(
    x = "Gene copy concentrations (GC mL<sup>-1</sup>)",
    y = NULL
  ) +
  theme(
    strip.text = element_markdown(face = "bold"),
    axis.title.x = element_markdown(),
    plot.margin = margin(l = 5, r = 15))

# ggsave("figures/qpcr-boxplots.png",
#        qpcr_boxplots,
#        dev = "png", dpi = 300,
#        width = 16, height = 8)


# Spearman correlation plot -----------------------------------------------

# Named vectors for plotting 
# Params to plot
these_params <- c("mcye_ndaf", "total_cyano", "ratio")

# Tidy labels
facet_labs <- list(
  analyte = c(
    mcye_ndaf = "*mcyE* (GC mL<sup>-1</sup>)",
    total_cyano = "Total cyanobacteria (GC mL<sup>-1</sup>)",
    ratio = "*mcyE*:Total cyanobacteria"
))

# Swap names and values for fct_recode()
param_labels <- setNames(names(facet_labs$analyte), facet_labs$analyte)
# Add a line break 
names(param_labels) <- names(param_labels) |> 
  str_replace( "(\\(.*\\))", "<br>\\1")

# Arrange data
spearman_data <- complete_data |> 
  mutate(
    # encode <lod as 0.01 to allow plotting without using another aesthetic
    total_cyano = if_else(total_cyano < lod_total_cyano, 0, log10(total_cyano)),
    mcye_ndaf = if_else(mcye_ndaf < lod_mcye_ndaf, 0, log10(mcye_ndaf)),
    ratio = mcye_ndaf / total_cyano
    )

# Locations with most MC-LR detections and max grab concentrations --------

most_detections_df <- spearman_data |> 
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

spearman_correlations <- spearman_data |>
  filter(!is.na(passive_mclr)) |> 
  select(location, all_of(these_params), contains("passive_mclr")) |> 
  # indicator column for MC-LR detections
  mutate(detected = passive_mclr > passive_mclr_lod) |> 
  select(-contains("mclr")) |> 
  group_by(location) |> 
  summarise(
    across(c(mcye_ndaf, total_cyano, ratio), ~ median(.x, na.rm = TRUE)),
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
                    c("mcye_ndaf", "total_cyano", "ratio")))

# Faceted jitter plot -----------------------------------------------------

# generate plot
plot_data <- spearman_data |>
  select(location, year,  contains("mcye"), contains("cyano"), ratio, contains("passive_mclr")) |> 
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
    !is.na(detected) 
  ) |> 
  left_join(most_detections_df, by = "location") |>
  # left_join(max_concentrations_df, by = "location") |> 
  mutate(
    # order the locations by the ordering found above
    label = factor(label, levels = rev(most_detections)),
    # label = factor(label, levels = rev(max_concentrations)),
    name = factor(name, levels = 
                    c("mcye_ndaf", "total_cyano", "ratio")),
    n_detections = factor(n_detections)
    # max = factor(max)
  ) 


spearman_plot_1 <- plot_data |> 
  filter(name == "mcye_ndaf") |> 
  ggplot(aes(label, value)) +
  facet_wrap(vars(name), nrow = 1,
             scales = "free_x",
             labeller = labeller(name = facet_labs$analyte)) +
  geom_jitter(aes(color = detected,  alpha = detected, size = detected), 
              show.legend = c(color = TRUE, alpha = FALSE, size = FALSE),
              width = 0.1) +
  coord_flip() +
  scale_color_manual(values = c("grey", detection_color)) +
  scale_alpha_manual(values = c(0.5, 0.9)) +
  scale_size_manual(values = c(2,3)) +
  scale_y_continuous(
    limits = NULL,
    breaks = c( 0, 2, 4, 6),
    labels =c("<LOD", "10<sup>2</sup>", "10<sup>4</sup>", "10<sup>6</sup>")) +
  guides(
    color = guide_legend(
      title = "MC-LR",
      position = "bottom",
      override.aes = list(size = 3)
    )
  ) +
  geom_richtext(
    data = spearman_correlations |> filter(name == "mcye_ndaf"),
    aes(x = Inf, y = Inf,
        label = paste("*&rho;* = ",spearman, sep = ""),
        hjust = 1, vjust =1,
        fontface = "bold")
  ) +
  theme(
    axis.text.y = element_markdown(hjust = 1),
    axis.text.x = element_markdown(),
    strip.text.x = element_markdown()
  ) +
  labs(
    x = NULL,
    y = NULL
  )

spearman_plot_2 <- plot_data |>
  filter(name == "total_cyano") |> 
  ggplot(aes(label, value)) +
  facet_wrap(vars(name), nrow = 1,
             scales = "free_x",
             labeller = labeller(name = facet_labs$analyte)) +
  geom_jitter(aes(color = detected,  alpha = detected, size = detected), 
              show.legend = c(color = TRUE, alpha = FALSE, size = FALSE),
              width = 0.1) +
  coord_flip() +
  scale_color_manual(values = c("grey", detection_color)) +
  scale_alpha_manual(values = c(0.5, 0.9)) +
  scale_size_manual(values = c(2,3)) +
  scale_y_continuous(
    limits = NULL,
    breaks = c( 0, 2, 4, 6),
    labels =c("<LOD", "10<sup>2</sup>", "10<sup>4</sup>", "10<sup>6</sup>")) +
  guides(
    color = guide_legend(
      title = "MC-LR",
      position = "bottom",
      override.aes = list(size = 3)
    )
  ) +
  geom_richtext(
    data = spearman_correlations |> filter(name == "total_cyano"),
    aes(x = Inf, y = Inf,
        label = paste("*&rho;* = ",spearman, sep = ""),
        hjust = 1, vjust =1,
        fontface = "bold")
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_markdown(),
    strip.text.x = element_markdown()
  ) +
  labs(
    x = NULL,
    y = NULL
  )

spearman_plot_3 <- plot_data |>
  filter(name == "ratio") |> 
  ggplot(aes(label, value)) +
  facet_wrap(vars(name), nrow = 1,
             scales = "free_x",
             labeller = labeller(name = facet_labs$analyte)) +
  geom_jitter(aes(color = detected,  alpha = detected, size = detected), 
              show.legend = c(color = TRUE, alpha = FALSE, size = FALSE),
              width = 0.1) +
  coord_flip() +
  scale_color_manual(values = c("grey", detection_color)) +
  scale_alpha_manual(values = c(0.5, 0.9)) +
  scale_size_manual(values = c(2,3)) +
  scale_y_continuous(
    limits = NULL,
    breaks = c( 0, 0.5, 1),
    labels =c("<LOD", "0.5", "1")) +
  guides(
    color = guide_legend(
      title = "MC-LR",
      position = "bottom",
      override.aes = list(size = 3)
    )
  ) +
  geom_richtext(
    data = spearman_correlations |> filter(name == "ratio"),
    aes(x = Inf, y = Inf,
        label = paste("*&rho;* = ",spearman, sep = ""),
        hjust = 1, vjust =1,
        fontface = "bold")
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_markdown(),
    strip.text.x = element_markdown()
  ) +
  labs(
    x = NULL,
    y = NULL
  )

spearman_plot <- spearman_plot_1 + spearman_plot_2 + spearman_plot_3 +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    strip.background = element_rect(color = "black", linewidth = 0.5)
    )
    
ggsave("figures/si-figures/qpcr-spearman-plot.png", spearman_plot, width = 14, height = 8)

