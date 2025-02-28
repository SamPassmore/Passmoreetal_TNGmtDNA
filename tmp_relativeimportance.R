library(phangorn)
library(phytools)
library(traitheritage)

mcct_tree = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.mcct.trees") 
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

trait = metadata$TNG
names(trait) = metadata$PhyID
distance = readRDS("processed_data/hikingtime_h.RDS")

keeps = metadata[metadata$ID %in% colnames(distance), c("PhyID", "ID")]
d2 = distance[rownames(distance) %in% metadata$ID, rownames(distance) %in% metadata$ID]

dimnames(d2) = list(keeps$PhyID[match(colnames(d2), keeps$ID)],
                    keeps$PhyID[match(colnames(d2), keeps$ID)])

trait = trait[names(trait) %in% rownames(d2)]

mcct_tree2 = keep.tip(mcct_tree, names(trait))

length(trait)
dim(d2)
length(mcct_tree2$tip.label)

test = traitheritage::trait_coevolution(mcct_tree2, trait, d2, 25, 12)

td = test$time_df
td$generation = as.numeric(gsub("g_", replacement = "", td$generation))
td$language_only = rowSums(td[,c("lang_dist", "lang_nodist")]) / td$denominator_sum

plot(x = td$generation,
     y = td$p_lang_dist, ylim = c(0,1),
     type = "l")
lines(x = td$generation, y = td$p_lang_nodist, col = "red")
lines(x = td$generation, y = td$p_nolang_dist, col = "green")
lines(x = td$generation, y = td$p_nolang_nodist, col = "blue")
lines(x = td$generation, y = td$language_only, col = "purple")

ans = ace(trait, mcct_tree)

plot(mcct_tree, type = "c", FALSE, label.offset = 1)
co <- c("blue", "yellow")
tiplabels(pch = 22, bg = co[as.numeric(trait + 1)], cex = 2, adj = 1)
nodelabels(pch = 19, col = rgb(0, 0, ans$ace), cex = 0.75)
