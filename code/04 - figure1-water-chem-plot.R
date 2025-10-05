
library(tidyverse)
library(dunn.test)
library(ggtext)

theme_set(theme_bw(base_size = 16))

detection_color <- "firebrick"
grab_col <- "#51c016"
passive_col <- "#4a0c6b"
year_colors <- c("#E69F00","#56B4E9", "#009E73")


complete_data <- read_csv(here::here("data-cleaned/data-complete.csv")) |> 
  mutate(
    year = as.factor(year),
    chloride = chloride / 1000,
    tn = tn * 1000) |> 
  rename(
    chloride_mgl = chloride,
    tn_ugl = tn,
    p_ugl = p,
    toc_mgl = toc,
    doc_mgl = doc)

# Kruskall Wallis and Dunn test -------------------------------------------

these_params <- c("p_ugl", "chloride_mgl", "tn_ugl", "colour", "doc_mgl", "toc_mgl")

# Kruskall Wallis
kw_results <- complete_data |> 
  mutate(location = factor(location)) |> 
  select(location, sample_date, all_of(these_params)) |> 
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

# Post-hoc dunn test with Benjamini Hochberg correction for multiple comparisons
dunn_results <- complete_data |> 
  mutate(location = factor(location)) |> 
  select(location, sample_date, all_of(these_params)) |> 
  pivot_longer(-c(location, sample_date)) |> 
  group_by(name) |> 
  nest() |> 
  mutate(
    dunn = map(data, ~ dunn.test(.x$value, .x$location, method = "bh") |>  as_tibble())) |> 
  ungroup() |> 
  select(name, dunn) |> 
  unnest(dunn) |> 
  arrange(name, P.adjusted)


write_csv(kw_results, "results/kw-test-water-chem.csv")
write_csv(dunn_results, "results/dunn_results-water-chem.csv")  


# Plot dunn test results --------------------------------------------------

facet_labs <- list(
  analyte = c(chloride_mgl = "Chloride (mg L<sup>-1</sup>)",
              p_ugl = "Total Phosphorus (&mu;g L<sup>-1</sup>)",
              tn_ugl = "Total Nitrogen (&mu;g L<sup>-1</sup>)",
              colour = "True Colour (ptCo)",
              phosphate_ugl = "Phosphate (&mu;g L<sup>-1</sup> as P)",
              toc_mgl = "Total Organic Carbon (mg L<sup>-1</sup>)",
              doc_mgl = "Dissolved Organic Carbon (mg L<sup>-1</sup>)"))

lake_labels <- c(
  "Shubenacadie Canal" = "SC",
  "Tower Road Reservoir" = "TRR",
  "Turtle Creek Reservoir" = "TCR"
)

# Get levels to allow filtering out redundant comparisons
location_levels <- sort(unique(complete_data$location))

# Wrangle
df <- dunn_results %>%
  #split into two columns
  separate(comparisons, into = c("lake1", "lake2"), sep = " - ") |> 
  mutate(
    lake1 = factor(lake1, levels = location_levels),
    lake2 = factor(lake2, levels = location_levels)
  ) |> 
  # remove redundant comparisons
  filter(as.integer(lake1) < as.integer(lake2))  |> 
  filter(P.adjusted < 0.05) |> #only keep significant comparison
  # rescale to make the aesthetics interpretable
  mutate(
    neglog_p = -log10(P.adjusted),
    alpha_val = scales::rescale(-log10(P.adjusted), to = c(0.4, 1)),  # Min alpha = 0.4
    size_val = scales::rescale(-log10(P.adjusted), to = c(4, 8)),       # Bubble size scale
    lake1 = recode(lake1, !!!lake_labels),
    lake2 = recode(lake2, !!!lake_labels)
    )


# Plot
dunn_pvalue_plot <- df |> 
  ggplot(aes(x = lake2, y = lake1)) +
  facet_wrap(vars(name),
             labeller = labeller(name = facet_labs$analyte)) +
  geom_point(aes(size = size_val, col = neglog_p)) +

  scale_color_viridis_c(
    option = "D",
    name = expression(-log[10](Adj.~italic(p))),
    limits = c(0, max(df$neglog_p, na.rm = TRUE))
  ) +
  scale_size_identity() +
  scale_alpha_identity() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    # panel.grid = element_blank(),
    strip.text = element_markdown(face = "bold"),
    axis.title = element_blank())
  
ggsave("figures/si-figures/dunn-pvalues.png", dunn_pvalue_plot, width = 17, height = 8)

# Plot lake water chemistry -----------------------------------------------

nutrient_plot <- complete_data |>
  select(location, year, all_of(these_params)) |> 
  pivot_longer(-c(location, year)) |> 
  filter(
    case_when(
      name == "p_ugl" & value > 250 ~ FALSE,
      name == "tn_ugl" & value > 1000 ~ FALSE,
      # name == "phosphate_ugl" & value > 50 ~ FALSE,
      TRUE ~ TRUE
      )
    ) |> 
  mutate(
    location =  factor(location, levels = location_levels),
    location = recode(location, !!!lake_labels),
    name = factor(name, levels = 
                    c("chloride_mgl", "p_ugl", "phosphate_ugl", 
                      "tn_ugl", "colour", "toc_mgl", "doc_mgl"))
  ) |> 
  ggplot(aes(location, value)) +
  facet_wrap(vars(name), scales = "free_y",
             labeller = labeller(name = facet_labs$analyte)) +
  geom_jitter(aes(color = year), width = 0.1, alpha = 0.4, size = 3) +
  geom_boxplot(alpha = 0) +
  expand_limits(y = 0) +
  geom_richtext(
    data = kw_results,
    aes(x = Inf, y = Inf,
        label = paste("Kruskal-Wallis *p* ", p_value_tidy, sep = ""),
        hjust = 1, vjust =1)
  ) +
  scale_color_manual(values = year_colors) +
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
    x = NULL,
    y = NULL
  ) +
  theme(
    strip.text = element_markdown(face = "bold"),
    axis.text.x = element_markdown(angle = 45, vjust = 1, hjust = 1),
    plot.margin = margin(l = 5, r = 15))

ggsave("figures/nutrients.png", nutrient_plot, width = 16, height = 8)

# Create tables of data ---------------------------------------------------

lake_chem_table_yearly <- complete_data |> 
  select(location, year, total_cyano, chloride_mgl, ph_field, p_ugl, tn_ugl, colour, toc_mgl) |> 
  pivot_longer(-c(location, year)) |> 
  group_by(location, year, name) |> 
  summarise(
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  ) |> 
  ungroup() |> 
  pivot_wider(
    names_from = name,
    values_from = c(median, sd)
  ) |> 
  mutate(across(where(is.numeric), ~ round(.x, 1)))

write_csv(lake_chem_table_yearly, "tables/lake_chem_table_yearly.csv")

lake_chem_table <- complete_data |> 
  select(location, total_cyano, chloride_mgl, ph_field, p_ugl, tn_ugl, colour, toc_mgl) |> 
  pivot_longer(-c(location)) |> 
  group_by(location, name) |> 
  summarise(
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  ) |> 
  ungroup() |> 
  pivot_wider(
    names_from = name,
    values_from = c(median, sd)
  )|> 
  mutate(across(where(is.numeric), ~ round(.x, 1)))

write_csv(lake_chem_table, "tables/lake_chem_table.csv")


# Water quality correlations ----------------------------------------------

# approach taken from: https://www.bryanshalloway.com/2020/06/03/tidy-2-way-column-combinations/
df_lists <- complete_data |>
  select(temp_field:tn_ugl, -contains("rdl"), -contains("mclr"), mcye_ndaf, total_cyano) |>
  summarise(across(where(is.numeric), list)) |> 
  pivot_longer(
    everything(),
    names_to = "var",
    values_to = "vector") |>
  mutate(
    var = fct_recode(var,
                     "Temp (C&deg;)" = "temp_field",
                     "DO (mg L<sup>-1</sup>)" = "do_field",
                     "Conductivity (&mu;S cm<sup>-1</sup>)" = "cond_field",
                     "pH" = "ph_field",
                     "TDS (mg L<sup>-1</sup>)" = "tds_field",
                     "*mcyE* (GC mL<sup>-1</sup>)" = "mcye_ndaf",
                     "Total cyanobacteria (GC mL<sup>-1</sup>)" = "total_cyano",
                     "Al (&mu;g L<sup>-1</sup>)" = "al",
                     "P (&mu;g L<sup>-1</sup>)" = "p_ugl",
                     "Fe (&mu;g L<sup>-1</sup>)" = "fe",
                     "Cl (mg L<sup>-1</sup>)" = "chloride_mgl",
                     "NO<sub>3</sub> (mg L<sup>-1</sup>)" = "nitrate",
                     "SO<sub>4</sub> (mg L<sup>-1</sup>)" = "sulfate",
                     "UV<sub>254</sub> (cm<sup>-1</sup>)" = "uv_254",
                     "Turbidity (NTU)" = "turbidity",
                     "True colour (PtCo)" = "colour",
                     "DOC (mg L<sup>-1</sup>)" = "doc_mgl",
                     "TOC (mg L<sup>-1</sup>)" = "toc_mgl",
                     "TN (&mu;g L<sup>-1</sup>)" = "tn_ugl"))
  
df_lists_comb <- expand(df_lists,
       nesting(var, vector),
       nesting(var2 = var, vector2 = vector)) |> 
  filter(var != var2) |> 
  arrange(var, var2) |> 
  mutate(vars = paste0(var, ".", var2)) |> 
  select(contains("var"), everything()) 

# Function to collapse redundant comparisons 
c_sort_collapse <- function(...){
  c(...) |> 
    sort() |> 
    str_c(collapse = ".")
}

pairwise_df <- df_lists_comb |> 
  mutate(vars = map2_chr(.x = var,
                         .y = var2,
                         .f = c_sort_collapse)) |> 
  distinct(vars, .keep_all = TRUE)

# Get pearson correlations and wrangle
correlations <- pairwise_df |> 
  mutate(
    cor_results = map2(vector, vector2, cor.test),
    correlation = map_dbl(cor_results, "estimate"),
    p_value = map_dbl(cor_results, "p.value")
  ) |> 
  select(var, var2, correlation, p_value) |> 
  arrange(desc(abs(correlation))) 

correlations |> 
  ggplot(aes(var, var2, size = abs(correlation), alpha = abs(correlation), color = correlation)) +
  geom_point() +
  scale_alpha_continuous(range = c(0.5, 0.95)) +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_distiller(palette = "PRGn", limits = c(-1, 1)) +
  guides(
    col = guide_colorbar(
      title = "Pearson <br> correlation"
      ),
    alpha = "none",
    size = "none"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.title = element_markdown(),
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    axis.text.y = element_markdown()
  )

ggsave("figures/si-figures/correlations.png",
       device = "png",
       dpi = 600,
       height = 6, width = 8)


