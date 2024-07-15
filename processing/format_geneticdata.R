# Cleaning genetic data
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(insight)
})

# Read in Raw data
data_file = 'data/PNG_PPP_Language GPS 270723.xlsx'

#### Format Sheet 1 ####
genetic_sheet1 = read_xlsx(data_file, sheet = 1, na = "ND")

# remove the last two rows being read in as they are just notes
genetic_sheet1 = genetic_sheet1[-c(696:697),]

genetic_sheet1 = genetic_sheet1 %>% 
  dplyr::select(NCBI_EHA_ID = "ID - NCBI/EGA" ,
                reference = "Reference",
                reference_sequence = "Reference Sequence ID",
                gender = "Gender",
                country = "Country",
                province = "Ind Province",
                district = "Ind District",
                llg = "Ind LLG", # LLG is a level of local government
                ward = "Ind Ward",
                village = "Ind Village/region**", 
                comments = "Comment **",
                region = "Population group / region",
                altitude = "Ind Altitude",
                latitude_ofsampling = "Ind Latitude*",
                longitude_ofsampling = "Ind Longitude",
                interview_language = "Ind Lang from interview",
                language_fromcoords = "Ind lang from gps",
                haplogroup = "Haplogroup",
                mutations_fromRCRS = "Mutations relative to Rcrs"
  )

# Remove completely empty rows
genetic_sheet1 = genetic_sheet1[rowSums(is.na(genetic_sheet1)) != ncol(genetic_sheet1),]

## Correct IDs
## Not all rows have IDs in this sheet
## Create unique IDs from the given ID and reference sequence ID
genetic_sheet1$ID = genetic_sheet1$reference_sequence
idx = is.na(genetic_sheet1$ID)
genetic_sheet1$ID[idx] = genetic_sheet1$NCBI_EHA_ID[idx]
genetic_sheet1$ID[genetic_sheet1$ID == "EGAS00001001247"] = genetic_sheet1$NCBI_EHA_ID[genetic_sheet1$ID == "EGAS00001001247"]

# If there is no language recorded, use the inferred language data 
cat("We infer the language for", 
    sum(is.na(genetic_sheet1$interview_language) & !is.na(genetic_sheet1$language_fromcoords)),
    "individuals in genetic sheet 1, from a total of",
    nrow(genetic_sheet1), "individuals")

# genetic_sheet1$interview_language[is.na(genetic_sheet1$interview_language)] = 
#   genetic_sheet1$language_fromcoords[is.na(genetic_sheet1$interview_language)]

#### Format Sheet 2 ####
genetic_sheet2 = read_xlsx(data_file, sheet = 2, skip = 3, na = "ND")
genetic_sheet2 = genetic_sheet2 %>% 
  dplyr::select(NCBI_EHA_ID = "ID - NCBI/EGA",
                reference = "Reference",
                reference_sequence = "Reference Sequence ID",
                gender = "Gender",
                country = "Country",
                province = "Ind Province",
                district = "Ind District",
                ward = "Ind Ward",
                village = "Ind Village/region**", 
                comments = "Comment **",
                region = "Population group / region",
                altitude = "Ind Altitude",
                latitude_ofsampling = "Ind Latitude*",
                longitude_ofsampling = "Ind Longitude",
                interview_language = "Ind Lang from Interview",
                language_fromcoords = "Ind lang from gps",
                haplogroup = "Haplogroup",
                mutations_fromRCRS = "Mutations relative to Rcrs",
                father_province = "Father Province",
                father_district = "Father District",
                father_llg = "Father LLG",
                father_ward = "Father Ward",
                father_village = "Father Village/region**",
                father_altitude = "Father Altitude",
                father_latitude = "Lat_Father",
                father_longitude = "Long_Father",
                mother_province = "Mother Province",
                mother_district = "Mother District",
                mother_llg = "Mother LLG",
                mother_ward = "Mother Ward",
                mother_village = "Mother Village/region**",
                mother_altitude = "Mother Altitude",
                mother_latitude = "Lat_Mother",
                mother_longitude = "Long_Mother",
                paternal_language = "Paternal Language Interview",
                maternal_language = "Maternal Language Interview",
                paternal_language_inferred = "Paternal Language inferred gps"	,
                maternal_language_inferred = "Maternal Language inferred gps")

## Not all rows have IDs in this sheet
## Create unique IDs from the reference sequence ID
genetic_sheet2$ID = genetic_sheet2$reference_sequence

#### Format Sheet 3 ####
genetic_sheet3 = read_xlsx(data_file, sheet = 3, skip = 1, na = "NA")

# remove the first row because it was only used for aesthetics
genetic_sheet3 = genetic_sheet3[-1,]

genetic_sheet3 = genetic_sheet3 %>%
  dplyr::select(ID = "Num_Form",
    interview_location = "Place_Sampling", # check these values are provinces
    gender = "Gender",
    village = "Place_Birth",
    district = "District_Birth",
    province = "Province_Birth",
    longitude_ofsampling = "Longitude",
    latitude_ofsampling = "Latitude",
    altitude = "Elevation",
    father_village = "F1_Place",
    father_district = "F1_District",
    father_province = "F1_Province",
    father_tribe = "F1_Tribe",
    father_clan = "F1_Clan",
    father_longitude = "F1_Longitude",
    father_latitude = "F1_Latitude",
    father_altitude = "F1_Elevation",
    mother_village = "M2_Place",
    mother_district = "M2_District",
    mother_province = "M2_Province",
    mother_tribe = "M2_Tribe",
    mother_clan = "M2_Clan",
    mother_longitude = "M2_Longitude",
    mother_latitude = "M2_Latitude",
    mother_altitude = "M2_Elevation",
    interview_language = "Language ind",
    interview_languagefamily = "Language_Family",
    father_language = "F1_Language_Glottolog",
    father_languagefamily = "F1_LanguageFamily_Glottolog",
    father_language_glottocode = "F1_LanguageCode_Glottolog",
    mother_language = "M2_Language_Glottolog",
    mother_language_glottocode = "M2_LanguageCode_Glottolog",
    mother_languagefamily = "M2_LanguageFamily_Glottolog",
    haplogroup = "Haplogroup",
    mutations_fromRCRS = "Mutations relative to Rcrs")

genetic_sheet3$interview_language[is.na(genetic_sheet3$interview_language)] = ""

genetic_sheet3$interview_language[genetic_sheet3$interview_language == "Kuman"] = "Kuman [kue]"
genetic_sheet3$interview_language[genetic_sheet3$interview_language == "Kuman and Pidgen"] = "Kuman [kue]"
genetic_sheet3$interview_language[genetic_sheet3$interview_language == "Kuman tokpies"] = "Kuman [kue]"
genetic_sheet3$interview_language[genetic_sheet3$interview_language == "Kuman/Pidgin/English"] = "Kuman [kue]"

### Format Sheet 4 ####
## The fourth sheet we received from Francois is the Daru data. We are ignoring data from the Daru sampling.

#### Testing data assumptions ####
## For each sheet ensure there are unique ids
n_distinct(genetic_sheet1$ID) == nrow(genetic_sheet1)
n_distinct(genetic_sheet2$ID) == nrow(genetic_sheet2)
n_distinct(genetic_sheet3$ID) == nrow(genetic_sheet3)

#### Bind rows ####
genetic_data = bind_rows(genetic_sheet1, genetic_sheet2, genetic_sheet3)
# missing language data should be an NA
genetic_data$interview_language[genetic_data$interview_language == ""] = NA

#### Picking the right language ####
## There are some genetic duplicates, of the duplicate we need to choose the individual with the most reliable metadata
## That is, we should aim to not use the inferred language data
dups = sum(duplicated(genetic_data$mutations_fromRCRS))

genetic_data$mutations_sets = as.numeric(factor(genetic_data$mutations_fromRCRS))

split_byset = split(genetic_data[,c("ID", "mutations_sets", "interview_language", "latitude_hybrid", "longitude_hybrid")], genetic_data$mutations_sets)
idx = sapply(split_byset, nrow) > 1
duplicates = do.call(rbind, split_byset[idx])

write.csv(duplicates, "processed_data/genetic_duplicates.csv")

dup_sets = sum(sapply(split_byset, nrow) > 1)

cat("There are", dups, "individuals with duplicate mtDNA mutations that group in to", 
    dup_sets, "distinct sets\n")

mismatches = sapply(split_byset, function(x) length(unique(x$interview_language)))

cat("Of the", dup_sets, "sets of duplicates", sum(mismatches > 1), "have mismatched language groupings\n")

for(i in unique(genetic_data$mutations_sets)){
  idx = which(genetic_data$mutations_sets == i)
  set = genetic_data[idx,]
  
  # If no languages have linguistic data, then use the inferred data
  if(all(is.na(set$interview_language))){
    set$interview_language = set$language_fromcoords
    # If the inferred languages are not all the same, then pick one at random
    if(length(unique(set$interview_language)) == 1 & nrow(set) > 1){
      drop = sample(idx, length(idx) - 1, replace = FALSE)
      genetic_data = genetic_data[-drop,]
    }
    next
  }

  # If all duplicate genomes speak the same langauge, do nothing
  if(length(unique(set$interview_language)) == 1){
    next
  }
  
  # If all languages are different, pick one at random
  if(all(!is.na(set$interview_language)) & !all(set$interview_language[1] == set$interview_language)){
    drop = sample(idx, length(idx) - 1, replace = FALSE)
    genetic_data = genetic_data[-drop,]
    print("These people have identical genetic records but different languages")
    print(set[,c("ID", "interview_language")])
    next
  }
  
  # If there are duplicate genomes but only one unique language, drop the others
  if(any(!is.na(set$interview_language))){
    drop = idx[is.na(set$interview_language)]
    genetic_data = genetic_data[-drop,]
    next
  }
}

# If there is no language recorded, use the inferred language data 
cat("We infer the language for", 
    sum(is.na(genetic_data$interview_language) & !is.na(genetic_data$language_fromcoords)),
    "individuals, from a total of",
    nrow(genetic_data), "individuals")

# genetic_sheet1$interview_language[is.na(genetic_sheet1$interview_language)] = 
#   genetic_sheet1$language_fromcoords[is.na(genetic_sheet1$interview_language)]


# Separately I collected the GPS coordinates of all Villages listed in the dataset. I merge that data in here
village_coords = read.csv('data/metadata_villagedata.csv')
village_coords = village_coords[,c("ID", "village", "LATITUDE", "LONGITUDE")]

genetic_data = left_join(genetic_data, village_coords, by = c("ID", "village"))

## In the metadata received from Brucato & Ricaut, there is no metadata for B00IEVO-PNG82 
## But some metadata exists for this individual in Pedro 2020 metadata. Why has this been dropped? 

## Make some assumptions based on varying questions
## If paternal / maternal language is missing, use the inferred language
### As it stands this is only one language. 
genetic_data$paternal_language[is.na(genetic_data$paternal_language)] = 
  genetic_data$paternal_language_inferred[is.na(genetic_data$paternal_language)]

genetic_data$maternal_language[is.na(genetic_data$maternal_language)] = 
  genetic_data$maternal_language_inferred[is.na(genetic_data$maternal_language)]

## Paternal and Maternal language can be inferred to be Mother and Father language
genetic_data$father_language[is.na(genetic_data$father_language)] = 
  genetic_data$paternal_language[is.na(genetic_data$father_language)]

genetic_data$mother_language[is.na(genetic_data$mother_language)] = 
  genetic_data$maternal_language[is.na(genetic_data$mother_language)]


#### More Tests ####
# Make M and N haplotype identifier
## Label M and N Haplogroups
## According to: http://phylotree.org/tree/index.htm
m_haplogroup = c("M", "G", "D", "Q", "C", "Z", "E")
n_haplogroup = c("N", "A", "X", "R", "B", "U", "O", "S", "I", "W", "Y", "P", "HV", "H", "V", "J", "T", "F", "K")
genetic_data$haplogroup_MN = ifelse(substr(genetic_data$haplogroup, 1, 1) %in% m_haplogroup, "M",
                           ifelse(substr(genetic_data$haplogroup, 1, 1) %in% n_haplogroup, "N",
                                  NA))

# are all haplotypes either M or N?
## we know that ID 525 does not have a haplotype in our current dataset
all(!is.na(genetic_data$haplogroup_MN[-890]))

cat("There are", nrow(genetic_data), "samples\n")
cat("There are", sum(genetic_data$haplogroup_MN == "M", na.rm = T), "M haplotype samples\n")
cat("There are", sum(genetic_data$haplogroup_MN == "N", na.rm = T), "N haplotype samples\n")

# Language data
cat("There are", n_distinct(genetic_data$interview_language), "languages in the sample\n")

# There are two individuals with the ID PM18 who have been sampled from different locations. 
# It is unclear why this difference exists, so we exclude both of them for now. 
## PM18 from Oro Province is in the Pedro 2020 Supplementary material 
genetic_data = genetic_data %>% 
  filter(ID != "PM18")

# all rows have IDS?
assertthat::assert_that(sum(is.na(genetic_data$ID)) == 0)
assertthat::assert_that(n_distinct(genetic_data$ID) == nrow(genetic_data))

# Put the ID columns first
genetic_data = relocate(genetic_data, "ID", 1)

#### Manual Fixes ####
## Manual fixes to some linguistic data
# For recoding purposes, it is easier to have missing data listed descriptively
genetic_data$interview_language[is.na(genetic_data$interview_language)] = "None given"
genetic_data$interview_language[genetic_data$interview_language == ""] = "None given"
# Interviewees gave their language as Angan, which is a language family. The closest Angan language to where the interviews
# were conducted was Baruya. These individuals are assigned to Baruya
genetic_data$interview_language[genetic_data$interview_language == "Angan"] = "Baruya [byr]"
# Nasioi is a family and Simeku is a language in that family. These individuals are assigned to Simeku
genetic_data$interview_language[genetic_data$interview_language == "Nasioi, Simeku"] = "Simeku [smz]"
# Interviewees gave their language as Baining, which is a language family. The closest Baining language to where the interviews
# were conducted was Ura These individuals are assigned to Ura
genetic_data$interview_language[genetic_data$interview_language == "Baining"] = "Ura [uro]"
# Kup-Minj is a dialect of Wahgi which is gives Kumai as an alternative name. 
genetic_data$interview_language[genetic_data$interview_language == "Kumai ?"] = "Wahgi [wgi]"
# Melpa [med], Wiru [wiu] are both given as interview languages. They are also listed as the paternal and maternal language respectively.
# We largely only have one language per person, and will make this a constraint of this analysis. We assign Wiru to this language. 
genetic_data$interview_language[genetic_data$interview_language == "Melpa [med], Wiru [wiu]"] = "Wiru [wiu]"
# Tokano [zuh], Dano [aso] are both given as interview languages. They are also listed as the paternal and maternal language respectively.
# We can only have one language per person, and will make this a constraint of this analysis. We assign Dano to these people. 
genetic_data$interview_language[genetic_data$interview_language == "Tokano [zuh], Dano [aso]"] = "Dano [aso]"
# Tiang [tbj], Tigak [tgc] are both given as interview languages. They are also listed as the paternal and maternal language respectively.
# We can only have one language per person, and will make this a constraint of this analysis. We assign Tigak to these people. 
genetic_data$interview_language[genetic_data$interview_language == "Tiang [tbj], Tigak [tgc]"] = "Tigak [tgc]"
# Mutu [tuc], Mape [mlh] are both given as interview languages. They are also listed as the paternal and maternal languages.
# We can only have one language per person, and will make this a constraint of this analysis. We assign Mutu to these people. 
genetic_data$interview_language[genetic_data$interview_language == "Mutu [tuc], Mape [mlh]"] = "Mutu [tuc]"
# Misima-Paneati [mpx], Suau [swp] are both given as interview languages. They are also listed as the paternal and maternal language respectively.
# We can only have one language per person, and will make this a constraint of this analysis. We assign Suau to these people. 
genetic_data$interview_language[genetic_data$interview_language == "Misima-Paneati [mpx], Suau [swp]"] = "Suau [swp]"
# Lele [lle], Tok Pisin are both listed as interview languages. Since Tok Pisin is a lingua franca, this response is listed as Lele
genetic_data$interview_language[genetic_data$interview_language == "Lele [lle], Tok Pisin"] = "Lele [lle]"
# Bola [bnp], Tok Pisin are both listed as interview languages. Since Tok Pisin is a lingua franca, this response is listed as Bola
genetic_data$interview_language[genetic_data$interview_language == "Bola [bnp], Tok Pisin"] = "Bola [bnp]"
# 3 individuals have their language recorded as "Gimi (Eastern Highlands) [gim]" and 8 as "Gimi [gim]".
# These appear to be the same languages, so they are changed to be the same
genetic_data$interview_language[genetic_data$interview_language == "Gimi (Eastern Highlands) [gim]"] = "Gimi [gim]"

## Two langauges cannot be matched to a recorded language variety
# Eyagu (1 person)
# Imani (17 people)
## Imani appears to be a valley / place where people have come from. It is unclear what language they speak. 

## One language is listed as TO EXCLUDE. We exclude it here
genetic_data = genetic_data[!genetic_data$interview_language == "TO EXCLUDE",]

# If the response to the languages question is listed as None given - recode to NA
genetic_data$interview_language[genetic_data$interview_language == "None given"] = NA

# check it worked
!any(genetic_data$interview_language == "TO EXCLUDE", na.rm = T)
!any(genetic_data$interview_language == "None given", na.rm = T)

#### Standardization ####
genetic_data$gender = ifelse(genetic_data$gender == "m", "M",
                            ifelse(genetic_data$gender == "f", "F",  genetic_data$gender))

# Is the number of unique IDs equal to the number of rows? (i.e. are there no duplicates)
n_distinct(genetic_data$ID) == nrow(genetic_data)

## Add TNG Language branches or langauge family  ####
## I have created a manual pairing of languages to language families / branches. 
tng = read.csv('data/language_families.csv')

genetic_data = left_join(genetic_data, tng, by = "interview_language")

#### Hybrid variables ####
### These variables contain a mixutre of observed and inferred data. Details for their hybridisation are below

# The location of the sample is taken from a mixture of sources. 
# Some locations are taken as the location of sampling
# Some locations are taken as the home village of the sampler 
# We prefereably use the location of sampling. 
genetic_data$latitude_hybrid = coalesce(genetic_data$latitude_ofsampling, genetic_data$LATITUDE)
genetic_data$longitude_hybrid = coalesce(genetic_data$longitude_ofsampling, genetic_data$LONGITUDE)

##### Write the raw data ####
# Comma seperated
write.csv(genetic_data, 'processed_data/genetic_samples.csv', row.names = FALSE, na = "")

# tab seperate for TreeViewer
picks = c("ID", "Branch_Family", "interview_language")
tsv = genetic_data[, picks]
tsv$Branch_Family[is.na(tsv$Branch_Family)] = "Missing"


write.table(
  tsv,
  'processed_data/genetic_samples.tsv',
  row.names = FALSE,
  na = "",
  sep = "\t",
  quote = FALSE
)
