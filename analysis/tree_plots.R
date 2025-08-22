## Haplogroup table
suppressPackageStartupMessages({
  library(ape)
  library(forcats)
  library(dplyr)
  library(kableExtra)
  library(knitr)
  library(janitor)
})

source("analysis/helper.R")
options(knitr.kable.NA = '-')

# Load Metadata
TNG_metadata = read.csv("processed_data/TNGData.csv")

## Add in duplicate metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% TNG_metadata$PhyID,] # subset to the duplicates only

TNG_metadata = rbind(TNG_metadata, TNG_dups)

TNG_metadata = TNG_metadata[TNG_metadata$interview_language != "Tok Pisin [tpi]",] 
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a metadata mistake. 
TNG_metadata = TNG_metadata[TNG_metadata$Branch_Family != "Baining",]

## Get MCCT Tree, and tree depth
mcct_tree = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.mcct.duptips.trees") 
mcct_depth = mean(adephylo::distRoot(mcct_tree))

TNG_metadata = TNG_metadata[TNG_metadata$PhyID %in% mcct_tree$tip.label,]

# Results of TNG binary trait
tng_results = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_probability_summary.csv")
tng_permutations = read.csv("results/TNG_20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips_permutedsummary.csv")

tng_results = tng_results[tng_results$generation < mcct_depth,]
tng_permutations = tng_permutations[tng_permutations$generation < mcct_depth,]

tng_summary = summarizer(tng_results)
tng_summary_p = summarizer(tng_permutations, permuted_results = TRUE)

# write.csv(tng_summary, "processed_data/TNG_summary.csv")
# write.csv(tng_summary_p, "processed_data/TNG_summary_p.csv")

# tng_summary = read.csv("processed_data/TNG_summary.csv")
# tng_summary_p = read.csv("processed_data/TNG_summary_p.csv")

TNG_df = TNG_metadata %>% dplyr::filter(PhyID %in% mcct_tree$tip.label)

tng_tree = keep.tip(mcct_tree, TNG_df$PhyID[!is.na(TNG_df$TNG)])
tree_data = full_join(tng_tree, TNG_df, by = c(label = "PhyID"))

tree_line_plot(tree_data, tng_summary, tng_summary_p, variable = "TNG") & theme(legend.position = "none")


#when does the data probability seperate from the permutation probability
# true_perm = left_join(tng_summary, tng_summary_p, by = c("generation"),
#                       suffix = c(".true", ".perm"))
# idx = which(true_perm$CI_low.true > true_perm$CI_high.perm)
# head(true_perm[idx,])
# 
# # # when does the data probability rejoin from the permutation probability
# idx = (true_perm$CI_low.true > true_perm$CI_high.perm)
# true_perm %>% filter(idx) %>% select(generation, CI_low.true, CI_high.perm) %>% View()
# 
# sum(true_perm$generation[idx] > 4225 & true_perm$generation[idx] < 10350)
# length(seq(from = 4225, to = 10350, by = 25))
