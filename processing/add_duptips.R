## This script adds duplicate tips back to the phylogeny

library(ape)
library(phangorn)

duplicates = read.csv("processed_data/genetic_duplicates.csv")
metadata = read.csv("processed_data//TNGData.csv")
trees  = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20241216_CBESbTest_LanguageSample/20241112_CBESbTest_LanguageSample-all_data.noburnin.trees") 
mcct.tree = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20241216_CBESbTest_LanguageSample/20241112_CBESbTest_LanguageSample-all_data.mcct.trees") 

tt = trees[[1]]

dup_sets = split(duplicates, duplicates$mutations_sets)

t2 = lapply(trees[1:10000], function(tt) {
  for (i in seq_along(dup_sets)) {
    set = dup_sets[[i]]

    if (any(set$PhyID %in% tt$tip.label)) {
      to_add = set[!set$PhyID %in% tt$tip.label, ]
      # should always be one
      exists = set[set$PhyID %in% tt$tip.label, ]
      
      tt = add.tips(tt, to_add$PhyID, exists$PhyID)
    }
  }
  tt
})


## MCCT Dup Tips
for (i in seq_along(dup_sets)) {
  set = dup_sets[[i]]
  
  if (any(set$PhyID %in% mcct.tree$tip.label)) {
    to_add = set[!set$PhyID %in% mcct.tree$tip.label, ]
    # should always be one
    exists = set[set$PhyID %in% mcct.tree$tip.label, ]
    
    mcct.tree = add.tips(mcct.tree, to_add$PhyID, exists$PhyID)
  }
}


write.nexus(t2, file = "/Volumes/PassmoreSSD/PapuanPast/20241216_CBESbTest_LanguageSample/20241112_CBESbTest_LanguageSample-all_data.noburnin.duptips.trees")
write.nexus(mcct.tree, file = "/Volumes/PassmoreSSD/PapuanPast/20241216_CBESbTest_LanguageSample/20241112_CBESbTest_LanguageSample-all_data.mcct.duptips.trees")

# Checks
# length(t2[[1]]$tip.label)
# # length(unique(duplicates$PhyID)) == length(duplicates$PhyID)
# all(duplicates$PhyID %in% t2[[1]]$tip.label)
# 
# which(!duplicates$PhyID %in% t2[[1]]$tip.label)
# 
# duplicates[!duplicates$PhyID %in% t2[[1]]$tip.label,]
# 
# sum(duplicates$PhyID %in% t2[[1]]$tip.label) / nrow(duplicates)
