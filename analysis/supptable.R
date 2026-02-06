library(ape)

## Get Metadata
TNG_metadata = read.csv("processed_data/TNGData.csv")
## Add in duplicate metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% TNG_metadata$PhyID,] # subset to the duplicates only
metadata = rbind(TNG_metadata, TNG_dups)
metadata = metadata[metadata$interview_language != "Tok Pisin [tpi]",]
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a metadata mistake. 
metadata = metadata[metadata$Branch_Family != "Baining",]


mcct_tree = read.nexus("/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.mcct.duptips.trees") 

md = metadata[metadata$PhyID %in% mcct_tree$tip.label,]

## Samples Table

supp = data.frame(ID = md$ID,
           Reference = md$citation,
           Language = md$interview_language,
           `Family/Branch` = md$Branch_Family,
           Haplogroup = md$haplogroup,
           mutations_fromRCRS = md$mutations_fromRCRS
)

supp$Reference[supp$Reference == "New"] = "André et al. (2021)"

supp = supp[order(supp$Reference),]

dim(supp)
sum(supp$Language == "Tok Pisin [tpi]")

write.csv(supp, "figures/supplementarytable_1.csv", row.names = FALSE, fileEncoding = "UTF-8")

tmp = which(md$interview_language == md$father_language | 
      md$interview_language == md$mother_language & 
      md$TNG == 1)

table(md$Branch_Family[tmp])

