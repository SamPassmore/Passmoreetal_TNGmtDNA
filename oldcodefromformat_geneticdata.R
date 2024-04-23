## Format languages codes
language_clean = genetic_data$interview_language %>% 
  str_remove_all("\\[.*\\]") %>%
  str_remove_all("\\?") %>%
  str_trim()

language_names = language_clean %>% 
  strsplit(., ", ") %>% 
  plyr::ldply(., rbind) 
language_names = apply(language_names, 2, str_trim)
colnames(language_names) = c("first_language", "second_language")
language_names = data.frame(language_names)
language_names$original_name = assigned_language

language_names = distinct(language_names)
language_names$iso639 = str_match(language_names$original_name, "\\[([a-z]{3})\\]")[,2]

# Pair language data with genetic data
# genetic_data = left_join(genetic_data, language_names, by = c("assigned_language" = "original_name")) %>% 
#   left_join(., languages, by = c("iso639" = "ISO639.single"), suffix = c('', ".language"))

# cat("We have the genetic samples of", nrow(genetic_data), "individuals")
# cat("We have a language paired with", sum(!is.na(genetic_data$iso639)), "individuals")

# remove rows with no genetic data
genetic_data = genetic_data[!is.na(genetic_data$mutations_fromRCRS),]

## Make Counts per location of genetic data
genetic_agg = genetic_data %>% 
  group_by(latitude_ofsampling, longitude_ofsampling, assigned_language) %>% 
  summarise(sample_count = n()) %>% 
  na.omit()

#### To Fix later ####
# There are two individuals with the ID PM18 who have been sampled from different locations. 
# It is unclear why this difference exists, so we exlcude both of them for now. 
genetic_data = genetic_data %>% 
  filter(ID != "PM18")

## Manual fixes to some linguistic data
# Interviewees gave their language as Angan, which is a language family. The closest Angan language to where the interviews
# were conducted was Baruya. These individuals are assigned to Baruya
genetic_data$assigned_language[genetic_data$assigned_language == "Angan"] = "Baruya [byr]"
# Nasioi is a family and Simeku is a language in that family. These individuals are assigned to Simeku
genetic_data$assigned_language[genetic_data$assigned_language == "Nasioi, Simeku"] = "Simeku [smz]"
# Interviewees gave their language as Baining, which is a language family. The closest Baining language to where the interviews
# were conducted was Ura These individuals are assigned to Ura
genetic_data$assigned_language[genetic_data$assigned_language == "Baining"] = "Ura [uro]"
# Kup-Minj is a dialect of Wahgi which is gives Kumai as an alternative name. 
genetic_data$assigned_language[genetic_data$assigned_language == "Kumai ?"] = "Wahgi [wgi]"
# Melpa [med], Wiru [wiu] are both given as interview languages. They are also listed as the paternal and maternal language respectively.
# We largely only have one language per person, and will make this a constraint of this analysis. We assign Wiru to this language. 
genetic_data$assigned_language[genetic_data$assigned_language == "Melpa [med], Wiru [wiu]"] = "Wiru [wiu]"
# Tokano [zuh], Dano [aso] are both given as interview languages. They are also listed as the paternal and maternal language respectively.
# We largely only have one language per person, and will make this a constraint of this analysis. We assign Dano to these people. 
genetic_data$assigned_language[genetic_data$assigned_language == "Tokano [zuh], Dano [aso]"] = "Dano [aso]"
# Tiang [tbj], Tigak [tgc] are both given as interview languages. They are also listed as the paternal and maternal language respectively.
# We largely only have one language per person, and will make this a constraint of this analysis. We assign Tigak to these people. 
genetic_data$assigned_language[genetic_data$assigned_language == "Tiang [tbj], Tigak [tgc]"] = "Tigak [tgc]"
# Mutu [tuc], Mape [mlh] are both given as interview languages. They are also listed as the paternal and maternal languages.
# We largely only have one language per person, and will make this a constraint of this analysis. We assign Mutu to these people. 
genetic_data$assigned_language[genetic_data$assigned_language == "Mutu [tuc], Mape [mlh]"] = "Mutu [tuc]"
# Misima-Paneati [mpx], Suau [swp] are both given as interview languages. They are also listed as the paternal and maternal language respectively.
# We largely only have one language per person, and will make this a constraint of this analysis. We assign Suau to these people. 
genetic_data$assigned_language[genetic_data$assigned_language == "Misima-Paneati [mpx], Suau [swp]"] = "Suau [swp]"
# Lele [lle], Tok Pisin are both listed as interview languages. Since Tok Pisin is a lingua franca, this response is listed as Lele
genetic_data$assigned_language[genetic_data$assigned_language == "Lele [lle], Tok Pisin"] = "Lele [lle]"
# Bola [bnp], Tok Pisin are both listed as interview languages. Since Tok Pisin is a lingua franca, this response is listed as Bola
genetic_data$assigned_language[genetic_data$assigned_language == "Bola [bnp], Tok Pisin"] = "Bola [bnp]"

## genetic sheet 4 ####
genetic_sheet4 = read_xlsx(data_file, sheet = 4, skip = 1, na = c("ND", "NA"))

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
                interview_langaugefamily = "Language_Family",
                father_language = "F1_Language_Glottolog",
                father_languagefamily = "F1_LanguageFamily_Glottolog",
                father_language_glottocode = "F1_LanguageCode_Glottolog",
                mother_language = "M2_Language_Glottolog",
                mother_language_glottocode = "M2_LanguageCode_Glottolog",
                mother_languagefamily = "M2_LanguageFamily_Glottolog",
                haplogroup = "Haplogroup",
                mutations_fromRCRS = "Mutations relative to Rcrs")

genetic_sheet4$ID = as.character(genetic_sheet4$ID)

## Fix some language labelling
genetic_sheet4$interview_language[genetic_sheet4$interview_language == "Anta"] = "Anta [tci]"
genetic_sheet4$interview_language[genetic_sheet4$interview_language == "KIwai"] = "Kiwaian [kjd]"
genetic_sheet4$interview_language[genetic_sheet4$interview_language == "Kiwai"] = "Kiwaian [kjd]"
