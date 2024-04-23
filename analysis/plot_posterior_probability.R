## Posterior Probability
library(dplyr)
library(tidyr)
library(bayestestR)
library(ggplot2)

posterior_languagefamily = readRDS("results/interview_languagefamily_mtDNA_N-papua_n_daru_posterior.RDS")
posterior_languagefamily = posterior_languagefamily[[1]]
names(posterior_languagefamily) = 1:length(posterior_languagefamily)

posterior_languagefamily = bind_rows(posterior_languagefamily, .id = "posterior")

posterior_languagefamily$generation2 = as.numeric(gsub("g_", replacement = "", posterior_languagefamily$generation))

## Summarise posterior
posterior_summary = posterior_languagefamily %>%
  filter(!is.na(clade_probability)) %>% # remove NA clades. Usually occur because of singleton clades
  group_by(generation2) %>%
  summarise(mean.prob = mean(clade_probability, na.rm = TRUE),
            ci_hdi = bayestestR::ci(clade_probability), # 95% credible interval
            n_matches = sum(numerator_sum)) %>% 
  unpack(ci_hdi)

# bayestestR::ci(posterior_languagefamily$clade_probability[posterior_languagefamily$generation2 == 50])

ggplot() + 
  geom_ribbon(data = posterior_summary, aes(y = mean.prob, x = generation2, ymin = CI_low, ymax = CI_high, group = 1),
              fill = "lightblue", alpha = 0.5) + 
  theme_classic() + 
  ylab("Probability") + xlab("Time (years)")
  



