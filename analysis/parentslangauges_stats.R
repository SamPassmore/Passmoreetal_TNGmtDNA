## Parents Language clades

metadata = read.csv("processed_data/TNGData.csv")
metadata = metadata[metadata$interview_language != "Tok Pisin [tpi]",]
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a metadata mistake. 
metadata = metadata[metadata$Branch_Family != "Baining",]
## Add in duplicate metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% metadata$PhyID,]
metadata = rbind(metadata, TNG_dups)
metadata = metadata[!is.na(metadata$PhyID),]

# Percent of parents who speak the same language
sum(metadata$father_language == metadata$mother_language, na.rm = TRUE) / 
  sum(complete.cases(metadata[,c("father_language", "mother_language")]))
sum(complete.cases(metadata[,c("father_language", "mother_language")]))

# Percent of parents who speak languages from the same language clade
sum(metadata$Branch_Family.father == metadata$Branch_Family.mother, na.rm = TRUE) / 
  sum(complete.cases(metadata[,c("Branch_Family.father", "Branch_Family.mother")]))
sum(complete.cases(metadata[,c("Branch_Family.father", "Branch_Family.mother")]))

# Percent of parents who speak languages from the same language clade
tng = metadata$TNG == 1

sum(metadata$TNG.father == metadata$TNG.mother, na.rm = TRUE) / 
  sum(complete.cases(metadata[,c("TNG.father", "TNG.mother")]))
sum(complete.cases(metadata[,c("TNG.father", "TNG.mother")]))

sum(metadata$TNG.father[tng] == metadata$TNG.mother[tng], na.rm = TRUE) / 
  sum(complete.cases(metadata[tng,c("TNG.father", "TNG.mother")]))




