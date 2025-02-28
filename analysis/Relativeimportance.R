suppressPackageStartupMessages({
  library(phangorn)
  library(phytools)
  library(traitheritage)
  library(gtools)
  library(dplyr)
})

n_trees = 50
n_perm = 5
trait_name = "TNG"

trees = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips.trees") 
trees = trees[(length(trees)-n_trees):length(trees)]
mcct_tree = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.mcct.duptips.trees") 
mcct_depth = mean(adephylo::distRoot(mcct_tree))

metadata = read.csv("processed_data/TNGData.csv")
metadata = metadata[metadata$interview_language != "Tok Pisin [tpi]",]
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a metadata mistake. 
metadata = metadata[metadata$Branch_Family != "Baining",]
## Add in duplicate metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% metadata$PhyID,]
metadata = rbind(metadata, TNG_dups)
metadata = metadata[!is.na(metadata$PhyID),]

metadata = metadata[metadata$PhyID %in% mcct_tree$tip.label,]
mcct_tree = keep.tip(mcct_tree, metadata$PhyID)

trait = metadata[,trait_name]
names(trait) = metadata$PhyID
distance = readRDS("processed_data/hikingtime_h.RDS")

# Ensure the distance matrix, the trait data, and the phylogeny match up.
keeps = metadata[metadata$ID %in% colnames(distance), c("PhyID", "ID")]
d2 = distance[rownames(distance) %in% metadata$ID, rownames(distance) %in% metadata$ID]

dimnames(d2) = list(keeps$PhyID[match(colnames(d2), keeps$ID)],
                    keeps$PhyID[match(colnames(d2), keeps$ID)])

trait = trait[names(trait) %in% rownames(d2)]

mcct_tree2 = keep.tip(mcct_tree, names(trait))
# prune the posterior
trees_p = lapply(trees, function(t) keep.tip(t, names(trait)))

# Checks
#all(names(trait) %in% rownames(d2))
#sapply(trees_p, function(t) all(t$tip_labels %in% names(trait)))


## Posterior Results
base_name = paste0(trait_name, "_", "20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips") 

# trait_coevo12 =
#   lapply(trees_p, function(t) {
#     traitheritage::trait_coevolution(t, trait, d2, 25, 12)
#   })
# saveRDS(trait_coevo12, paste0("results/", base_name, "_12hours.RDS"))
# 
# trait_coevo24 =
#   lapply(trees_p, function(t) {
#     traitheritage::trait_coevolution(t, trait, d2, 25, 24)
#   })
# saveRDS(trait_coevo24, paste0("results/", base_name, "_24hours.RDS"))
# 
# trait_coevo48 =
#   lapply(trees_p, function(t) {
#     traitheritage::trait_coevolution(t, trait, d2, 25, 48)
#   })
# saveRDS(trait_coevo48, paste0("results/", base_name, "_48hours.RDS"))
# 
# ## Summarised Results
# results_summary = function(df){
#   ld_summary = lapply(df, function(x) x$time_df)
#   ld_summary = bind_rows(ld_summary, .id = 'posterior')
#   ld_summary %>% 
#     group_by(generation) %>% 
#     summarise(p_lang_dist_av = mean(p_lang_dist, na.rm = TRUE),
#               p_lang_nodist_av = mean(p_lang_nodist, na.rm = TRUE),
#               p_nolang_dist_av = mean(p_nolang_dist, na.rm = TRUE),
#               p_nolang_nodist_av = mean(p_nolang_nodist, na.rm = TRUE)) %>% 
#     mutate(generation =  as.numeric(gsub("g_", replacement = "", generation))) %>% 
#     reshape2::melt(id = c("generation"))
# }
# 
# ld12 = results_summary(trait_coevo12)
# ld24 = results_summary(trait_coevo24)
# ld48 = results_summary(trait_coevo48)
# 
# write.csv(ld12, file = paste0("results/", trait_name, "_", base_name, "_12summary.csv"))
# write.csv(ld24, file = paste0("results/", trait_name, "_", base_name, "_24summary.csv"))
# write.csv(ld48, file = paste0("results/", trait_name, "_", base_name, "_48summary.csv"))
# 
# rm(trait_coevo12)
# rm(trait_coevo24)
# rm(trait_coevo48)

cat("Starting Permutations...\n")
## Permuted Results
set.seed(28367)

cat("12 hours...\n")
trait_coevo12_p =
  lapply(1:n_perm, function(p) {
    lapply(trees_p, function(t) {
      tp = gtools::permute(trait)
      traitheritage::trait_coevolution(t, tp, d2, 25, 12)$time_df
    })
  }
)
trait_coevo12_p = lapply(trait_coevo12_p, function(x) {
  names(x) = NULL
  bind_rows(x, .id = "permutation")
  })
trait_coevo12_p = bind_rows(trait_coevo12_p, .id = "posterior")
write.csv(trait_coevo12_p, paste0("results/", base_name, "_perm_12hours.csv"), row.names = FALSE)

cat("24 hours...\n")
set.seed(28367)
trait_coevo24_p =
  lapply(1:n_perm, function(p) {
    lapply(trees_p, function(t) {
      tp = gtools::permute(trait)
      traitheritage::trait_coevolution(t, tp, d2, 25, 24)$time_df
    })
  }
)
trait_coevo24_p = lapply(trait_coevo24_p, function(x) {
  names(x) = NULL
  bind_rows(x, .id = "permutation")
})
trait_coevo24_p = bind_rows(trait_coevo24_p, .id = "posterior")
write.csv(trait_coevo24_p, paste0("results/", base_name, "_perm_24hours.csv"), row.names = FALSE)

cat("48 hours...\n")
set.seed(28367)
trait_coevo48_p =
  lapply(1:n_perm, function(p) {
    lapply(trees_p, function(t) {
      tp = gtools::permute(trait)
      traitheritage::trait_coevolution(t, tp, d2, 25, 48)$time_df
    })
  }
)
trait_coevo48_p = lapply(trait_coevo48_p, function(x) {
  names(x) = NULL
  bind_rows(x, .id = "permutation")
})
trait_coevo48_p = bind_rows(trait_coevo48_p, .id = "posterior")
write.csv(trait_coevo48_p, paste0("results/", base_name, "_perm_48hours.csv"), row.names = FALSE)



