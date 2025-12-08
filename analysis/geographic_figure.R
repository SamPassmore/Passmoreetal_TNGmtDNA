library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#### Parameters ####
png_colors <- list(
  # Flag Colors (National Identity)
  flag = c("#D50032", "#FFCC00", "#000000", "#006400"),
  
  # Natural Colors (Landscapes & Biodiversity)
  nature = c("#007A33", "#005FAD", "#4A4A4A", "#FF6F3C"),
  
  # Cultural & Traditional Colors
  culture = c("#D27D2C", "#8B5A2B", "#40E0D0", "#E63946")
)

#### Functions ####
summarise_posterior = function(df){
  df %>% 
    select(generation, p_lang_nodist) %>% 
    filter(generation != 0) %>% 
    group_by(generation) %>% 
    pivot_longer(cols = c("p_lang_nodist"), names_to = "variable", 
                 values_to = "value") %>% 
    summarise(p_lang_nodist_probability = mean(value, na.rm = TRUE),
              ci_hdi = bayestestR::ci(value)) %>% 
    ungroup() %>% 
    tidyr::unpack(ci_hdi) 
}

#### Data ####
# 48 hours
ld48_df = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_coevoTNG_48h.csv")
ld48_df_p = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_coevoTNG_48h_permu.csv")

ld48_df_s1 = ld48_df %>% 
  select(generation, p_lang_dist, p_lang_nodist, p_nolang_dist, p_nolang_nodist) %>% 
  filter(generation != 0) %>% 
  group_by(generation) %>% 
  pivot_longer(cols = c("p_lang_dist", "p_lang_nodist", "p_nolang_dist", "p_nolang_nodist"),
               names_to = "variable", values_to = "value") %>% 
  group_by(variable, generation) %>%
  summarise(value = mean(value))

# 24 hours
ld24_df = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_coevoTNG_24h.csv")
ld24_df_p = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_coevoTNG_24h_permu.csv")

# 12 hours
ld12_df = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_coevoTNG_12h.csv")
ld12_df_p = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_coevoTNG_12h_permu.csv")


plot_time = function(df, perm, title) {
  ggplot() +
    geom_ribbon(
      data = df,
      aes(ymin = CI_low, ymax = CI_high, x = generation),
      fill = "#fae484",
      alpha = 0.5
    ) + 
    geom_line(
      data = df,
      aes(x = generation, y = p_lang_nodist_probability),
      col = "#FFCC00"
    ) +
    geom_ribbon(
      data = perm,
      aes(ymin = CI_low, ymax = CI_high, x = generation),
      fill = "lightblue",
      alpha = 0.5
    ) +
    geom_line(
      data = perm,
      aes(x = generation, y = p_lang_nodist_probability),
      col = "#00274D"
    ) +
    xlim(c(25, 25000)) +
    ylab("Probability of L & ¬D") +
    xlab("Time") +
    scale_x_continuous(
      breaks = c(5000, 10000, 15000, 20000),
      labels = c("5k", "10k", "15k", "20k"),
      limits = c(0, 25000)
    ) + 
    ggtitle(title) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text = element_text(size = 12))
}

ld48_df_s = summarise_posterior(ld48_df)
ld48_df_ps = summarise_posterior(ld48_df_p)

ld24_df_s = summarise_posterior(ld24_df)
ld24_df_ps = summarise_posterior(ld24_df_p)

ld12_df_s = summarise_posterior(ld12_df)
ld12_df_ps = summarise_posterior(ld12_df_p)

p1 = ggplot(ld48_df_s1, aes(x = generation, y = value, fill = variable), bg = "black") + 
  geom_area(position = "stack") +
  scale_fill_manual(name = "Types of Pairs",
                    labels = c("L & D",
                               "L & ¬D",
                               "¬L & D",
                               "¬L & ¬D"),
                    values = c(png_colors$flag)) + 
  xlim(0, 25000) +
  scale_x_continuous(
    breaks = c(5000, 10000, 15000, 20000),
    labels = c("5k", "10k", "15k", "20k"),
    limits = c(0, 25000)
  ) + 
  ylab("Cumulative probability") + 
  ggtitle("48 Hours Walking Distance") + 
  theme_classic() + 
  theme(legend.position = "bottom",
        axis.title = element_blank())

p12 = plot_time(ld12_df_s, ld12_df_ps, "12 Hours")
p24 = plot_time(ld24_df_s, ld24_df_ps, "24 Hours")
p48 = plot_time(ld48_df_s, ld48_df_ps, "48 Hours")

p3 = p12 + p24 + p48 + plot_layout(axis_titles = "collect")

p1 / p3 + plot_layout(heights = c(4, 2))

ggsave(filename = "figures/geographic_figure.png", width = 210, height = 297/2, units = "mm")
