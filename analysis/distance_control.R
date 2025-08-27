## Controlling for distance
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source('analysis/helper.R')

summarise_distance = function(df){
  df = lapply(df, function(x) x$time_df)
  df = bind_rows(df, .id = "state")
  
  tmp2 = df %>% 
    mutate(generation =  as.numeric(gsub("g_", replacement = "", generation))) %>% 
    group_by(generation) %>%
    summarise(lang_dist = mean(p_lang_dist, na.rm = TRUE),
              lang_dist_ci_hdi = bayestestR::ci(p_lang_dist),
              lang_nodist = mean(p_lang_nodist, na.rm = TRUE),
              lang_nodist_ci_hdi = bayestestR::ci(p_lang_nodist),
              nolang_dist = mean(p_nolang_dist, na.rm = TRUE),
              nolang_dist_ci_hdi = bayestestR::ci(p_nolang_dist),
              nolang_nodist = mean(p_nolang_nodist, na.rm = TRUE),
              nolang_nodist_ci_hdi = bayestestR::ci(p_nolang_nodist)
    ) %>% 
    tidyr::unpack(c(lang_dist_ci_hdi, lang_nodist_ci_hdi, nolang_dist_ci_hdi, nolang_nodist_ci_hdi),
                  names_sep = "_")
  
  tidyr::pivot_longer(tmp2, cols = c("lang_dist", "lang_nodist", "nolang_dist", "nolang_nodist"), 
                      names_to = "state")
}

## Summarise Permutations
sum_perm = function(df){
  df %>% 
    mutate(generation =  as.numeric(gsub("g_", replacement = "", generation))) %>% 
    filter(generation < 50000) %>% 
    mutate(posterior = rep(1:51, each = (2000 * 5))) %>% 
    group_by(generation, posterior) %>%
    summarise(wntree_cladeprob = mean(p_lang_nodist, na.rm = TRUE)) %>% 
    ungroup() %>% 
    group_by(generation) %>% 
    summarise(p_lang_nodist_probability = mean(wntree_cladeprob, na.rm = TRUE),
              ci_hdi = bayestestR::ci(wntree_cladeprob)) %>% 
    ungroup() %>% 
    tidyr::unpack(ci_hdi)  
}


# Papua New Guinea Color Palette
png_colors <- list(
  # Flag Colors (National Identity)
  flag = c("#D50032", "#FFCC00", "#000000", "#006400"),
  
  # Natural Colors (Landscapes & Biodiversity)
  nature = c("#007A33", "#005FAD", "#4A4A4A", "#FF6F3C"),
  
  # Cultural & Traditional Colors
  culture = c("#D27D2C", "#8B5A2B", "#40E0D0", "#E63946")
)
                    

ld12_df = readRDS("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_12hours.RDS")
ld24_df = readRDS("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_24hours.RDS") 
ld48_df = readRDS("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_48hours.RDS") 

ld12_df_summary = summarise_distance(ld12_df)
ld24_df_summary = summarise_distance(ld24_df)
ld48_df_summary = summarise_distance(ld48_df)

ld12_df_p = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_perm_12hours.csv") 
ld24_df_p = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_perm_24hours.csv") 
ld48_df_p = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_perm_48hours.csv") 

ld12_df_ps = sum_perm(ld12_df_p)
ld24_df_ps = sum_perm(ld24_df_p)
ld48_df_ps = sum_perm(ld48_df_p)

## why am I getting the same output
ld12_df_p %>% group_by(generation) %>% summarise(mean(lang_nodist, na.rm = TRUE))
lapply(ld12_df, function(x) x$time_df) %>%  bind_rows(., .id = "state") %>%  group_by(generation) %>% summarise(mean(lang_nodist, na.rm = TRUE))



p1 = ggplot(ld48_df_summary, aes(x = generation, y = value, fill = state), bg = "black") + 
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
  theme_classic() + 
  theme(legend.position = "bottom",
        axis.title = element_blank())

plot_time = function(df, perm, title) {
  ggplot() +
    geom_line(
      data = df[df$state == "lang_nodist", ],
      aes(x = generation, y = value, fill = state),
      col = "#C8102E"
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

head(ld12_df_summary$value[ld12_df_summary$state == "lang_nodist"], 10)
head(ld12_df_p$p_lang_nodist, 10)

p12 = plot_time(ld12_df_summary, ld12_df_ps, "12 Hours")
p24 = plot_time(ld24_df_summary, ld24_df_ps, "24 Hours")
p48 = plot_time(ld48_df_summary, ld48_df_ps, "48 Hours")

p3 = p12 + p24 + p48 + plot_layout(axis_titles = "collect")

p1 / p3 + plot_layout(heights = c(4, 2))
