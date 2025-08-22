## By TNG Lanuguage clade 
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Summarise the data for Angan, Kainantu-Goroka, and Chimbu-Waghi and save it
## Because it takes some time to run.

# tng_branch = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_probability_bytrait.csv")
# tb_permuted = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_permutedbytrait.csv")
# 
# # get permuted results I need
# angan_dfp = summarizer(tb_permuted[tb_permuted$state == "Angan",], permuted_results = TRUE)
# kaingo_dfp = summarizer(tb_permuted[tb_permuted$state == "Kainantu-Goroka",], 
#                         permuted_results = TRUE)
# chimbu_dfp = summarizer(tb_permuted[tb_permuted$state == "Chimbu-Wahgi",],
#                         permuted_results = TRUE)

# # save results
# write.csv(angan_dfp, "processed_data/angan_permuted_results.csv")
# write.csv(kaingo_dfp, "processed_data/kainantugoroka_permuted_results.csv")
# write.csv(chimbu_dfp, "processed_data/chimbuwahgi_permuted_results.csv")

# read results
angan_dfp = read.csv("processed_data/angan_permuted_results.csv")
kaingo_dfp = read.csv("processed_data/kainantugoroka_permuted_results.csv")
chimbu_dfp = read.csv("processed_data/chimbuwahgi_permuted_results.csv")

# tng_branch = read.csv("results/Branch_Family_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_probability_bytrait.csv")
tng_branch = read.csv("results/TNG_Branch_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_probability_bytrait.csv")

tngbranch_summary = tng_branch %>%
  filter(!is.na(clade_probability)) %>% # remove NA clades. Usually occur because of singleton clades
  dplyr::filter(generation <= mcct_depth) %>%
  dplyr::filter(
    !state %in% c(
      "Ok-Oksapmin",
      "Madang",
      "Mailuan",
      "Eleman",
      "Indo-European",
      "Ndu",
      "Sepik",
      "South Bougainville",
      "Wiru"
    )
  ) %>%
  group_by(generation, state) %>%
  summarise(
    mean_clade_probability = mean(clade_probability, na.rm = TRUE),
    ci_hdi = bayestestR::ci(clade_probability, method = "HDI"),
    # 95% credible interval
    n_matches = sum(numerator_sum)
  ) %>%
  tidyr::unpack(ci_hdi)

# write.csv(tngbranch_summary, "processed_data/tng_bytraitresults2.csv")
# tngbranch_summary = read.csv("processed_data/tng_bytraitresults2.csv")

tngbranch_summary$state = factor(tngbranch_summary$state,
                                 levels = c("Arafundi", "Austronesian", 
                                            "Bosavi", "Eastern Trans Fly",
                                            "Lower Sepik",
                                            "Pahoturi", "Yam", "NotTNG", "Angan", 
                                            "Chimbu-Wahgi", "Enga-Kewa-Huli", 
                                            "Kainantu-Goroka", "Kiwaian",
                                            "Suki-Gogodala"))
# Reduce categories for plot
tngbranch_summary$state = fct_collapse(tngbranch_summary$state,
                                       `Non-TNG` = c("Arafundi", "Austronesian", 
                                                     "Bosavi", "Eastern Trans Fly",
                                                     "Lower Sepik",
                                                     "Pahoturi", "Yam"),
                                       Angan = c("Angan"),
                                       `Chimbu-Wahgi` = "Chimbu-Wahgi", 
                                       `Enga-Kewa-Huli` = "Enga-Kewa-Huli", 
                                       `Kainantu-Goroka` = "Kainantu-Goroka", 
                                       `Kiwaian` = "Kiwaian",
                                       `Suki-Gogodala` = "Suki-Gogodala")

tngbranch_summary = tngbranch_summary %>% group_by(state, generation) %>%
  summarise(mean_clade_probability = sum(mean_clade_probability), 
            CI_low = sum(CI_low),
            CI_high = sum(CI_high))

tng_cols = colorRampPalette(c("#FFFFFF", "#F4CCCC", "#E69999", "#D96666", "#CC3333", 
                              "#C8102E", "#A50D24", "#800A1A", "#550610"))

p1 = ggplot() + 
  geom_area(data = tngbranch_summary,
            aes(x = generation, y = mean_clade_probability, fill = state), colour = "black") + 
  scale_x_continuous(breaks = c(5000, 10000, 15000, 20000), 
                     labels = c("5k", "10k", "15k", "20k"), limits = c(0, 25000)) +
  scale_fill_manual(values = c("#FFE066", rev(tng_cols(6)))) +
  theme_classic() + 
  ylab("Probability of a Shared Language History") + 
  # ylim(c(0, 0.4)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", linewidth=0.5),
        legend.position = "bottom", 
        legend.title = element_blank(),
        axis.title = element_blank()) 

## Proportions of Angan and Kainantu-Goroka
angan_df = tngbranch_summary[tngbranch_summary$state == "Angan",]
angan_df$CI_high[angan_df$CI_high > 0.3] = 0.3 # clip the error bars to make the graph readable
kaingo_df = tngbranch_summary[tngbranch_summary$state == "Kainantu-Goroka",]
chimbu_df = tngbranch_summary[tngbranch_summary$state == "Chimbu-Wahgi",]
chimbu_df$CI_high[chimbu_df$CI_high > 0.3] = 0.3 # clip the error bars to make the graph readable

p_angan = ggplot() + 
  geom_ribbon(
    data = angan_df,
    aes(ymin = CI_low, ymax = CI_high, x = generation),
    fill = "#F4CCCC",
    alpha = 0.5
  ) + 
  geom_line(
    data = angan_df,
    aes(x = generation, y = mean_clade_probability),
    col = "#C8102E"
  ) +
  geom_ribbon(
    data = angan_dfp,
    aes(ymin = CI_low, ymax = CI_high, x = generation),
    fill = "lightblue",
    alpha = 0.5
  ) +
  geom_line(
    data = angan_dfp,
    aes(x = generation, y = mean_clade_probability),
    col = "#00274D"
  ) +
  xlim(c(25, 25000)) +
  ylim(c(0, 0.33)) +
  ylab("Probability of a\n Shared Language") +
  xlab("Time") +
  scale_x_continuous(
    breaks = c(5000, 10000, 15000, 20000),
    labels = c("5k", "10k", "15k", "20k"),
    limits = c(0, 25000)
  ) + 
  ggtitle("Angan") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  ),
  axis.text = element_text(size = 10))

p_kaingo = ggplot() + 
  geom_ribbon(
    data = kaingo_df,
    aes(ymin = CI_low, ymax = CI_high, x = generation),
    fill = "#F4CCCC",
    alpha = 0.5
  ) + 
  geom_line(
    data = kaingo_df,
    aes(x = generation, y = mean_clade_probability),
    col = "#C8102E"
  ) +
  geom_ribbon(
    data = kaingo_dfp,
    aes(ymin = CI_low, ymax = CI_high, x = generation),
    fill = "lightblue",
    alpha = 0.5
  ) +
  geom_line(
    data = kaingo_dfp,
    aes(x = generation, y = mean_clade_probability),
    col = "#00274D"
  ) +
  xlim(c(25, 25000)) +
  # ylim(c(0, 0.33)) +
  ylab("Probability of a\n Shared Language") +
  xlab("Time") +
  scale_x_continuous(
    breaks = c(5000, 10000, 15000, 20000),
    labels = c("5k", "10k", "15k", "20k"),
    limits = c(0, 25000)
  ) + 
  ggtitle("Kainantu-Goroka") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  ),
  axis.text = element_text(size = 10))


p_cw = ggplot() + 
  geom_ribbon(
    data = chimbu_df,
    aes(ymin = CI_low, ymax = CI_high, x = generation),
    fill = "#F4CCCC",
    alpha = 0.5
  ) + 
  geom_line(
    data = chimbu_df,
    aes(x = generation, y = mean_clade_probability),
    col = "#C8102E"
  ) +
  geom_ribbon(
    data = chimbu_dfp,
    aes(ymin = CI_low, ymax = CI_high, x = generation),
    fill = "lightblue",
    alpha = 0.5
  ) +
  geom_line(
    data = chimbu_dfp,
    aes(x = generation, y = mean_clade_probability),
    col = "#00274D"
  ) +
  xlim(c(25, 25000)) +
  ylim(c(0, 0.33)) +
  ylab("Probability of a\n Shared Language") +
  xlab("Time") +
  scale_x_continuous(
    breaks = c(5000, 10000, 15000, 20000),
    labels = c("5k", "10k", "15k", "20k"),
    limits = c(0, 25000)
  ) + 
  ggtitle("Chimbu-Wahgi") +
  theme_classic() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    hjust = 1
  ),
  axis.text = element_text(size = 10))

p_bottom = (p_angan + p_cw + p_kaingo) + plot_layout(axis_titles = "collect")
p_figure = p1 / p_bottom

ggsave(plot = p_figure, filename = "figures/by_clade.png", height = 5)

# dates
angan_timedepth = angan_df %>% 
  left_join(., angan_dfp, by = c("generation"), suffix = c(".true", ".perm")) %>% 
  mutate(diff = CI_low.true - CI_high.perm) %>% 
  filter(diff > 0.01) %>% 
  summarise(start_end = formatC(range(generation), format = "d", big.mark = ","))

kaingor_timedepth = kaingo_df %>% 
  left_join(., kaingo_dfp, by = c("generation"), suffix = c(".true", ".perm")) %>%
  mutate(diff = CI_low.true - CI_high.perm) %>% 
  filter(diff > 0.01) %>%
  summarise(start_end = formatC(range(generation), format = "d", big.mark = ","))

cw_timedepth = chimbu_df %>% 
  left_join(., chimbu_dfp, by = c("generation"), suffix = c(".true", ".perm")) %>%
  mutate(diff = CI_low.true - CI_high.perm) %>% 
  filter(diff > 0.01) %>%
  summarise(start_end = formatC(range(generation), format = "d", big.mark = ","))
