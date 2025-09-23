library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

#### Functions #### 
summarise_posterior = function(df){
  df %>% 
    select(generation, clade_probability) %>% 
    filter(generation != 0) %>% 
    group_by(generation) %>% 
    pivot_longer(cols = c("clade_probability"), names_to = "variable", 
                 values_to = "value") %>% 
    summarise(mean_clade_probability = mean(value, na.rm = TRUE),
              ci_hdi = bayestestR::ci(value, ci = 0.99, method = "HDI")) %>% 
    ungroup() %>% 
    tidyr::unpack(ci_hdi) 
}

plot_subgroup = function(df, df_p, title){
  ggplot() + 
    geom_ribbon(
      data = df,
      aes(ymin = CI_low, ymax = CI_high, x = generation),
      fill = "#F4CCCC",
      alpha = 0.5
    ) + 
    geom_line(
      data = df,
      aes(x = generation, y = mean_clade_probability),
      col = "#C8102E"
    ) +
    geom_ribbon(
      data = df_p,
      aes(ymin = CI_low, ymax = CI_high, x = generation),
      fill = "lightblue",
      alpha = 0.5
    ) +
    geom_line(
      data = df_p,
      aes(x = generation, y = mean_clade_probability),
      col = "#00274D"
    ) +
    ylab("Probability of a\n Shared Language") +
    xlab("Time") +
    scale_x_continuous(
      breaks = c(5000, 10000, 15000, 20000, 25000, 30000),
      labels = c("5k", "10k", "15k", "20k", "25k", "30k"),
      limits = c(0, 35000)
    ) +
    ylim(c(0, 1)) + 
    ggtitle(title) +
    theme_classic() +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text = element_text(size = 10))
}

#### Data ####
angan = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Angan_probability_bytrait.csv")
angan_p = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Angan_permutedsummary.csv")

chimbu = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Chimbu-Wahgi_probability_bytrait.csv")
chimbu_p = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Chimbu-Wahgi_permutedsummary.csv")

kainantu = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Kainantu-Goroka_probability_bytrait.csv")
kainantu_p = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Kainantu-Goroka_permutedsummary.csv")

suki = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Suki-Gogodala_probability_bytrait.csv")
suki_p = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Suki-Gogodala_permutedsummary.csv")

enga = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Enga-Kewa-Huli_probability_bytrait.csv")
enga_p = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Enga-Kewa-Huli_permutedsummary.csv")

kiwaian = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Kiwaian_probability_bytrait.csv")
kiwaian_p = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_Kiwaian_permutedsummary.csv")

#### Summarise #### 
angan_s = summarise_posterior(angan)
angan_ps = summarise_posterior(angan_p)

chimbu_s = summarise_posterior(chimbu)
chimbu_ps = summarise_posterior(chimbu_p)

kainantu_s = summarise_posterior(kainantu)
kainantu_ps = summarise_posterior(kainantu_p)

suki_s = summarise_posterior(suki)
suki_ps = summarise_posterior(suki_p)

enga_s = summarise_posterior(enga)
enga_ps = summarise_posterior(enga_p)

kiwaian_s = summarise_posterior(kiwaian)
kiwaian_ps = summarise_posterior(kiwaian_p)

#### Plot #### 
angan_plt = plot_subgroup(angan_s, angan_ps, "Angan (59)")
chimbu_plt = plot_subgroup(chimbu_s, chimbu_ps, "Chimbu-Wahgi (63)")
kainantu_plt = plot_subgroup(kainantu_s, kainantu_ps, "Kainantu-Goroka (n = 107)")
suki_plt = plot_subgroup(suki_s, suki_ps, "Suki-Gogodala (6)")
enga_plt = plot_subgroup(enga_s, enga_ps, "Enga-Kewa-Huli (21)")
kiwaian_plt = plot_subgroup(kiwaian_s, kiwaian_ps, "Kiwaian (5)")

top = (kainantu_plt + chimbu_plt)
bottom = (angan_plt + enga_plt)


(top + plot_layout(axes = "collect")) /
  (bottom + plot_layout(axes = "collect"))

ggsave(filename = "figures/language_subgroup.png", width = 210, height = 297/2, units = "mm")

## Dates of seperation
#when does the data probability seperate from the permutation probability
date_sep = function(df, perm){
  true_perm = left_join(df, perm, by = c("generation"),
                        suffix = c(".true", ".perm"))
  
  idx = which(true_perm$CI_low.true < true_perm$CI_high.perm) 
  true_perm[idx,] %>% 
    select(generation, CI_low.true, CI_high.perm) 
}

# date_sep(angan_s, angan_ps) %>% View() # 25875
# date_sep(chimbu_s, chimbu_ps) %>% View() # 14550
# date_sep(kainantu_s, kainantu_ps) %>% View() # 14050
# date_sep(suki_s, suki_ps) %>% View() # 5725
# date_sep(enga_s, enga_ps) %>% View() # 2050
# date_sep(kiwaian_s, kiwaian_ps) %>% View() # 0
# 
