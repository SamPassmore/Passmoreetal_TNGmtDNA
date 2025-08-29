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
                reference = "Reference Sequence ID",
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
                mutations_fromRCRS = "Mutations relative to Rcrs",
                citation = "Reference"
  )

# Remove completely empty rows
genetic_sheet1 = genetic_sheet1[rowSums(is.na(genetic_sheet1)) != ncol(genetic_sheet1),]

## Correct IDs
## Not all rows have IDs in this sheet
## Create unique IDs from the given ID and reference sequence ID
genetic_sheet1$ID = genetic_sheet1$reference
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
                maternal_language_inferred = "Maternal Language inferred gps",
                citation = "Reference")

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
    mutations_fromRCRS = "Mutations relative to Rcrs",
    citation = "Reference")

genetic_sheet3$interview_language[is.na(genetic_sheet3$interview_language)] = ""

genetic_sheet3$interview_language[genetic_sheet3$interview_language == "Kuman"] = "Kuman [kue]"
genetic_sheet3$interview_language[genetic_sheet3$interview_language == "Kuman and Pidgen"] = "Kuman [kue]"
genetic_sheet3$interview_language[genetic_sheet3$interview_language == "Kuman tokpies"] = "Kuman [kue]"
genetic_sheet3$interview_language[genetic_sheet3$interview_language == "Kuman/Pidgin/English"] = "Kuman [kue]"

### Format Sheet 4 ####
## The fourth sheet we received from Francois is the Daru data. We are ignoring data from the Daru sampling.
genetic_sheet4 = read_xlsx(data_file, sheet = 4, skip = 1, na = "NA")

# remove the first row because it was only used for aesthetics
genetic_sheet4 = genetic_sheet4[-1,]

genetic_sheet4 = genetic_sheet4 %>%
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
                mutations_fromRCRS = "Mutations relative to Rcrs",
                citation = "Reference")

genetic_sheet4$ID = as.character(genetic_sheet4$ID)

# Agob doesn't have an ISO code, so we use the Agob-Ende-Kawam ISO code
genetic_sheet4[genetic_sheet4 == "Agob"] = "Agob [kit]"
# Anta doesn't have an ISO code, so we use the Anta-Komna-Wara-Were-Kema ISO code
genetic_sheet4[genetic_sheet4 == "Anta"] = "Anta [tci]"
genetic_sheet4[genetic_sheet4 == "Arammba"] = "Arammba [stk]"
genetic_sheet4[genetic_sheet4 == "Ende"] = "Ende [end]"
# Metadata doesnt speciy what Kiwaian varity, but father speaks Southern Kiwai, so we label it this way
genetic_sheet4[genetic_sheet4 == "English/Kiwai"] = "Kiwai [kjd]"
# Gidra is an alternative name for Wipi on Glottolog
genetic_sheet4[genetic_sheet4 == "Gidra"] = "Wipi [gdr]"
genetic_sheet4[genetic_sheet4 == "Gizra"] = "Gizra [tof]"
genetic_sheet4[genetic_sheet4 == "Gododala"] = "Gogodala [ggw]"
genetic_sheet4[genetic_sheet4 == "Gogodala"] = "Gogodala [ggw]"
# one person lists gogodala and agob. We pick the first since we will only allow on language per person
genetic_sheet4[genetic_sheet4 == "Gogodala; Agob"] = "Gogodala [ggw]"
# ISO code is for Idi-Taemae
genetic_sheet4[genetic_sheet4 == "Idi"] = "Idi [idi]"
genetic_sheet4[genetic_sheet4 == "Idi (Yidi)"] = "Idi [idi]"
genetic_sheet4[genetic_sheet4 == "Idi / Taeme"] = "Idi [idi]"
genetic_sheet4[genetic_sheet4 == "Idi / Taeme"] = "Idi [idi]"
# Kawam doesn't have an ISO code, so we use the Agob-Ende-Kawam ISO code
genetic_sheet4[genetic_sheet4 == "Kawam (Dukop)"] = "Kawam [kit]"
# Not sure what this langauge is supposed to be
genetic_sheet4[genetic_sheet4 == "Kiunum"] = NA
## Assuming that Kiwai speakers are Southern Kiwai speakers
genetic_sheet4[genetic_sheet4 == "Kiwai"] = "Kiwai [kjd]"
genetic_sheet4[genetic_sheet4 == "KIwai"] = "Kiwai [kjd]"
genetic_sheet4[genetic_sheet4 == "Kiwai, Southern"] = "Kiwai [kjd]"
genetic_sheet4[genetic_sheet4 == "Kunja"] = ""
# Magiam is an alternative name for Kanch on Glottolog
genetic_sheet4[genetic_sheet4 == "Magiam"] = "Kánchá [pep]"
genetic_sheet4[genetic_sheet4 == "Makayam"] = "Makayam [aup]"
genetic_sheet4[genetic_sheet4 == "Nama"] = "Nama [nmx]"
genetic_sheet4[genetic_sheet4 == "Náma"] = "Nama [nmx]"
genetic_sheet4[genetic_sheet4 == "Nambo"] = "Nambo [ncm]"
genetic_sheet4[genetic_sheet4 == "Nen"] = "Nen [nqn]"
genetic_sheet4[genetic_sheet4 == "Nmbo"] = "Nambo [ncm]"
genetic_sheet4[genetic_sheet4 == "Suki"] = "Suki [sui]"
# Taeme doesnt have its own ISO code so we use Idi-Taeme
genetic_sheet4[genetic_sheet4 == "Taeme"] = "Taeme [idi]"
# Anta doesn't have an ISO code, so we use the Anta-Komna-Wara-Were-Kema ISO code
genetic_sheet4[genetic_sheet4 == "Wara"] = "Wára [tci]"
# Anta doesn't have an ISO code, so we use the Anta-Komna-Wara-Were-Kema ISO code
genetic_sheet4[genetic_sheet4 == "Wára"] = "Wára [tci]"
genetic_sheet4[genetic_sheet4 == "Wartha"] = "Warta Thuntai [gnt]"

#### Testing data assumptions ####
## For each sheet ensure there are unique ids
n_distinct(genetic_sheet1$ID) == nrow(genetic_sheet1)
n_distinct(genetic_sheet2$ID) == nrow(genetic_sheet2)
n_distinct(genetic_sheet3$ID) == nrow(genetic_sheet3)
n_distinct(genetic_sheet4$ID) == nrow(genetic_sheet4)

#### Bind rows ####
genetic_data = bind_rows(genetic_sheet1, genetic_sheet2, genetic_sheet3, genetic_sheet4)
# missing language data should be an NA
genetic_data$interview_language[genetic_data$interview_language == ""] = NA
genetic_data$mother_language[genetic_data$mother_language == ""] = NA
genetic_data$father_language[genetic_data$father_language == ""] = NA

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
## According to Gajdusek and Alpers (1972) Imani is a dialect of Baruya consisting of around 300 speakers (1969). We relabel them as Baruya
genetic_data$interview_language[genetic_data$interview_language == "Imani"] = "Baruya [byr]"
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

## One langauges cannot be matched to a recorded language variety
# Eyagu (1 person)

#### Father's language
## Simplify this note, but it is not clear what language Kora is. 
genetic_data$father_language[genetic_data$father_language == "Likely Kora of the Chimbu area."] = "Kora"
genetic_data$father_language[genetic_data$father_language == "Kuman"] = "Kuman [kue]"

#### Mother's language
genetic_data$mother_language[genetic_data$mother_language == "Asaro’o"] = "Asaro’o [mtv]"
# Bel appears to be a Austronesian clade. Class this as missing
genetic_data$mother_language[genetic_data$mother_language == "Bel"] = NA
# As with individual speakers, it is unclear what language Eyagu is.
genetic_data$mother_language[genetic_data$mother_language == "Eyagu ?"] = NA
genetic_data$mother_language[genetic_data$mother_language == "Gende"] = "Gende [gaf]"
## adding consistent labelling for Kuman
genetic_data$mother_language[genetic_data$mother_language == "Kuman"] = "Kuman [kue]"
## adding consistent labelling for North Wahgi
genetic_data$mother_language[genetic_data$mother_language == "Wahgi, North"] = "Wahgi, North [whg]"
genetic_data$mother_language[genetic_data$mother_language == "Kambai/Tapei [afp]"] = "Tapei [afp]"
genetic_data$mother_language[genetic_data$mother_language == "Biaru-Waria [WER]"] = "Biaru-Waria [wer]"

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
n_distinct(genetic_data$mutations_fromRCRS) == nrow(genetic_data)

## I have created a manual pairing of languages to language families / branches. 
tng = read.csv('data/language_families.csv')

## Adding Branch data
genetic_data = left_join(genetic_data, tng, by = "interview_language")
## Adding Branch data for father / mother
genetic_data = left_join(genetic_data, tng[,c("interview_language", "TNG", "Branch_Family", "TNG_Branch")], 
               by = c("father_language" = "interview_language"),
               suffix = c("", ".father"))

genetic_data = left_join(genetic_data, tng[,c("interview_language", "TNG", "Branch_Family", "TNG_Branch")], 
               by = c("mother_language" = "interview_language"),
               suffix = c("", ".mother"))

cat("After removing problematic or repeated samples there are", nrow(genetic_data), "mitogenomic samples.")

#### Hybrid variables ####
### These variables contain a mixture of observed and inferred data. Details for their hybridisation are below

# The location of the sample is taken from a mixture of sources. 
# Some locations are taken as the location of sampling
# Some locations are taken as the home village of the sampler 
# We prefereably use the location of sampling. 
genetic_data$latitude_hybrid = coalesce(genetic_data$latitude_ofsampling, genetic_data$LATITUDE)
genetic_data$longitude_hybrid = coalesce(genetic_data$longitude_ofsampling, genetic_data$LONGITUDE)

## Add PhyID
genetic_data = genetic_data %>% 
  mutate(PhyID = paste0(ID, "_", haplogroup),
       PhyID = gsub(" ", "", PhyID, fixed = TRUE))

## Merge and subset to TNG data
drop_provinces = c("West New Britain", 
                   "Bougainville",
                   "East New Britain",
                   "New Ireland")

cat("We drop", sum(genetic_data$province %in% drop_provinces) , "samples for being in provinces that are not of interest")
cat("We drop", sum(str_detect(genetic_data$haplogroup, "B4a")) , "samples for being of Austronesian origin")
cat("We drop", sum(is.na(genetic_data$interview_language)) , "samples for not having linguistic data")

# These individuals are identified as out of scope.
drop_ids = c("AY289086", "AY289092", "AY289087",
             "AY289084", "AY289088", "AY289091",
             "AY289081", "AY289090", "AY289085",
             "AY289082", "AY289079", "AY289078")

TNG = genetic_data %>% 
  filter(!str_detect(haplogroup, "B4a")) %>%  # remove Austronesian
  filter(!province %in% drop_provinces) %>% 
  filter(!ID %in% drop_ids) %>%  
  filter(!is.na(interview_language)) %>% nrow()

cat(sum(!is.na(TNG$interview_language)),
    "individuals provide language data, from a total of",
    nrow(genetic_data), "individuals")

# If there is no language recorded, use the inferred language data 
cat("We could infer the language for", 
    sum(is.na(TNG$interview_language) & !is.na(TNG$language_fromcoords)),
    "individuals, from a total of",
    nrow(TNG), "individuals. But we do not currently")

## Identify Austronesian languages
TNG$Austronesian = TNG$interview_language %in% c("Abadi [kbt]",
                                                 "Adzera [adz]",
                                                 "Kilivila [kij]",
                                                 "Maisin [mbq]",
                                                 "Mutu [tuc]",
                                                 "Nyindrou [lid]",
                                                 "Suau [swp]",
                                                 "Waima [rro]",
                                                 "Kuanua [ksd]")


#### Picking the right language ####
## There are some genetic duplicates, of the duplicate we need to choose the individual with the most reliable metadata
## That is, we should aim to not use the inferred language data
dups = sum(duplicated(TNG$mutations_fromRCRS))

TNG$mutations_sets = as.numeric(factor(TNG$mutations_fromRCRS))

split_byset = split(TNG[,c("PhyID", "mutations_sets", "interview_language", "latitude_hybrid", "longitude_hybrid")], TNG$mutations_sets)
idx = sapply(split_byset, nrow) > 1
duplicates = do.call(rbind, split_byset[idx])

duplicate_data = TNG[TNG$PhyID %in% duplicates$PhyID,]

write.csv(duplicates, "processed_data/genetic_duplicates.csv")
write.csv(duplicate_data, "processed_data/genetic_duplicate_data.csv", row.names = FALSE)

dup_sets = sum(sapply(split_byset, nrow) > 1)

cat("There are", dups, "individuals with duplicate mtDNA mutations that group in to",
    dup_sets, "distinct sets\n")

mismatches = sapply(split_byset, function(x) length(unique(x$interview_language)))

cat("Of the", dup_sets, "sets of duplicates", sum(mismatches > 1), "have mismatched language groupings\n")
record = list()
for(i in unique(TNG$mutations_sets)){
  idx = which(TNG$mutations_sets == i)
  set = TNG[idx,]
  
  # If no languages have linguistic data, then drop
  if(all(is.na(set$interview_language))){
    drop = idx
    TNG = TNG[-drop,]
    next
    # set$interview_language = set$language_fromcoords
    # If the inferred languages are not all the same, then pick one at random
    # if(length(unique(set$interview_language)) == 1 & nrow(set) > 1){
    #   drop = sample(idx, length(idx) - 1, replace = FALSE)
    #   genetic_data = genetic_data[-drop,]
    # }
  }
  
  # If all duplicate genomes speak the same language, then take the first one
  if(length(unique(set$interview_language)) == 1 & nrow(set) > 1){
    drop = idx[-1]
    TNG = TNG[-drop,]
    next
  }
  
  # If all languages are different, pick one at random
  if(any(is.na(set$interview_language)) | !all(set$interview_language[1] == set$interview_language)){
    drop = sample(idx, length(idx) - 1, replace = FALSE)
    TNG = TNG[-drop,]
    print("These people have identical genetic records but different languages")
    record = append(record, list(set[,c("ID", "interview_language")]))
    next
  }
  
  # If there are duplicate genomes but only one unique language, drop the others
  if(any(!is.na(set$interview_language)) & nrow(set) > 1){
    drop = idx[!is.na(set$interview_language)] # remove individuals with missing language data
    
    l_sample = ifelse(length(drop) == 1, 1, length(drop) - 1)
    drop = sample(drop, l_sample, replace = FALSE) # select langauge at random
    TNG = TNG[-drop,]
    next
  }
}

nrow(TNG)
n_distinct(TNG$mutations_fromRCRS) == nrow(TNG)

#### Write the raw data ####
# All Formatted data
write.csv(genetic_data, 'processed_data/genetic_samples.csv', row.names = FALSE, na = "")

# Only mainland PNG data
write.csv(TNG, "./processed_data/TNGData.csv", row.names = FALSE)
