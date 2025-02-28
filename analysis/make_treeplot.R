## make densitree

library(ggtree)
library(ape)
library(treeio)
library(ggplot2)
library(ggnewscale)


metadata = read.csv("processed_data/TNGData.csv")
metadata = metadata[metadata$interview_language != "Tok Pisin [tpi]",]
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a metadata mistake. 
metadata = metadata[metadata$Branch_Family != "Baining",]
## Add in duplicate metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% metadata$PhyID,]
metadata = rbind(metadata, TNG_dups)
metadata = metadata[!is.na(metadata$PhyID),]

posterior = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips.trees")
mcct = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.mcct.duptips.trees")

## First a densitree
p1 = ggdensitree(posterior[(length(posterior)-500):length(posterior)], alpha = 0.01, colour = "steelblue") 
ggsave(plot = p1, filename = "figures/20250106_LanguageSample_CBES_Relaxed_Constrained_densitree.png")

## Second the consensus tree w/ labels
first_tree = treeio::full_join(mcct, metadata, by = c(label = "PhyID"))

p2 = ggtree(first_tree, size = 0.1) + 
  geom_tippoint(aes(col = factor(TNG)), shape = "I", size = 4) +
  scale_x_continuous(breaks = c(10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000),
                     labels = c("10k", "20k", "30k", "40k", "50k", "60k", "70k", "80k")) +
  scale_fill_discrete(name = NULL) #+ 
  # layout_dendrogram() 
  # theme(legend.position = "bottom", 
  #       axis.ticks.x = element_line(),
  #       axis.line.x = element_line(),
  #       axis.text.x = element_text(),
  #       axis.text.y = element_blank(),
  #       panel.grid.major = element_line()) + 
  # theme_tree2()

p2

lang = as.data.frame(metadata[,c("PhyID", "TNG", "Branch_Family")])
rownames(lang) = metadata$PhyID
ggplot(lang, aes(x = TNG, y = PhyID)) + 
  geom_tile()
