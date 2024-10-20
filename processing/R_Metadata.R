## NOTE: All data formatting done in this script is now done in format_geneticdata.R

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

metadata = read.csv("processed_data/genetic_samples.csv")
load("figures/base_map.Rdata")  

duplicates = read.csv("processed_data/genetic_duplicates.csv")

# remove duplicate 
metadata2 = metadata %>% dplyr::filter(!ID %in% duplicates$ID)

## All TNG Speakers
## Remove Austronesian Haplogroups
## Remove Non-mainland provinces
drop_provinces = c("West New Britain", 
                   "Bougainville",
                   "East New Britain",
                   "New Ireland")

TNG_test = metadata %>% 
  filter(!str_detect(haplogroup, "B4a")) %>% # remove austronesians Austronesian
  filter(!province %in% drop_provinces) %>% 
  mutate(PhyID = paste0(ID, "_", haplogroup),
         PhyID = gsub(" ", "", PhyID, fixed = TRUE))

p = base_map + 
  geom_count(data = TNG_test, aes(x = longitude_hybrid, 
                                  y = latitude_hybrid),
             fill = "#FF5732",
             pch = 4, alpha = 0.5) + 
  ylab("Latitude") + xlab("Longitude") + 
  guides(col=guide_legend(title=element_blank())) +
  scale_color_hue(labels = c("Not Daru" = "Existing", "Daru" = "New")) + 
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank())

p

write.csv(TNG_test, "NCI/20240917_RAll/TNGData.csv")

## Read in Genetic data
gen_data = Biostrings::readDNAMultipleAlignment("processed_data/papua_completeset.fasta")
# subset to those of interest
gen_subset = gen_data@unmasked[TNG_test$PhyID]

# write.csv(TNG_test[,c("ID", "latitude_hybrid", "longitude_hybrid")], "processed_data/ArcGIS_LOCATIONS.csv",
#           row.names = FALSE)
Biostrings::writeXStringSet(gen_subset, filepath = 'NCI/20240910_CBESbModelTest_ALL/all.fasta')

Branch_test = metadata %>% 
  filter(TNG == 1)

p = base_map + 
  geom_count(data = Branch_test, aes(x = longitude_hybrid, 
                                  y = latitude_hybrid),
             fill = "#FF5732",
             pch = 4, alpha = 0.5) + 
  ylab("Latitude") + xlab("Longitude") + 
  guides(col=guide_legend(title=element_blank())) +
  scale_color_hue(labels = c("Not Daru" = "Existing", "Daru" = "New")) + 
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank())

p

# write.csv(Branch_test, "NCI/20240917_RAll/TNGBranchData.csv", row.names = FALSE)