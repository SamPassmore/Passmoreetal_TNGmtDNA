## This script adds duplicate tips back to the phylogeny

library(ape)
library(phangorn)

duplicates = read.csv("processed_data/genetic_duplicates.csv")

trees  = read.nexus("data/full_model-all.noburnin.trees")

tt = trees[[1]]

dup_sets = split(duplicates, duplicates$mutations_sets)

i = 2

t2 = lapply(trees[1:100], function(tt) {
  for (i in seq_along(dup_sets)) {
    set = dup_sets[[i]]

    if (any(set$PhyID %in% tt$tip.label)) {
      to_add = set[!set$PhyID %in% tt$tip.label, ]
      # should always be one
      exists = set[set$PhyID %in% tt$tip.label, ]
      
      print(exists$PhyID)
      print(to_add$PhyID)
      
      print(" ")
      tt = add.tips(tt, to_add$PhyID, exists$PhyID)
    }
  }
  tt
})


for (i in seq_along(dup_sets)) {
  set = dup_sets[[i]]
  
  print(sum(set$PhyID %in% trees[[1]]$tip.label))
} 
  
  
  
    
# }