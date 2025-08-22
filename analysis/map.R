suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
})

load('figures/base_map.Rdata')

# Load Metadata
TNG_metadata = read.csv("processed_data/TNGData.csv")

## Add in duplicate metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% TNG_metadata$PhyID,] # subset to the duplicates only

TNG_metadata = rbind(TNG_metadata, TNG_dups)

TNG_metadata = TNG_metadata[TNG_metadata$interview_language != "Tok Pisin [tpi]",] 
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a metadata mistake. 
TNG_metadata = TNG_metadata[TNG_metadata$Branch_Family != "Baining",]

# remove data not used in the phylogeny
mcct_tree = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.mcct.duptips.trees") 

TNG_metadata = TNG_metadata[TNG_metadata$PhyID %in% mcct_tree$tip.label,]

location_summary = TNG_metadata %>%
  mutate(long = round(longitude_hybrid, 1),
         lat =  round(latitude_hybrid, 1)) %>% 
  group_by(long, lat, TNG) %>% 
  summarise(n = n()) %>% 
  filter(long < 150) %>% 
  filter(!is.na(TNG))

p = base_map +
  geom_jitter(data = location_summary, aes(
    x = long,
    y = lat,
    size = n,
    col = factor(TNG),
    shape = factor(TNG)
  ),
  alpha = 0.5
  ) + 
  scale_size_continuous(breaks = c(10, 15, 20, 25)) + 
  scale_color_manual(labels = c("Not TNG", "TNG"), values = c("#FFCD00", "#C8102E")) + 
  scale_shape_manual(labels = c("Not TNG", "TNG"), values = c(17, 18)) + 
  guides(size = guide_legend(title ="N Samples", nrow = 1),
         col = guide_legend(title = element_blank()),
         shape = guide_legend(title = element_blank())) + 
  theme(legend.position = "bottom", 
        axis.title = element_blank())


ggsave(filename = "figures/map_figure.png")
