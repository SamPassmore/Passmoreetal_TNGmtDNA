suppressPackageStartupMessages({
  library(readxl)
  library(stringr)
  library(dplyr)
  library(assertthat)
})

glottolog = read.csv('https://raw.githubusercontent.com/glottolog/glottolog-cldf/master/cldf/languages.csv')

genetic_data = read.csv('processed_data/genetic_samples.csv', na.strings = "")

language_codes = str_extract_all(genetic_data$assigned_language, "(?<=\\[).+?(?=\\])") 
language_codes[sapply(language_codes, function(x) length(x) == 0)] = NA
language_codes = unlist(language_codes)

data = cbind(genetic_data, assigned_language_code = language_codes)

## Join glottolog
data = left_join(data, glottolog, by = c("assigned_language_code" = "ISO639P3code"), suffix = c("", ".glottolog"))

### Add language family data when missing ####
# Pawaia is a language with two dialects. The language family code is also the language code.
data$Family_ID[data$assigned_language == "Pawaia [pwa]"] = "pawa1255"
# Pele-Ata is a language with two dialects. The language family code is also the language code.
data$Family_ID[data$assigned_language == "Pele-Ata [ata]"] = "pele1245"


cat("There are", n_distinct(data$assigned_language_code) - 2, "languages identified in the sample.\n")
cat("There are", n_distinct(data$Family_ID) - 2, "Language families identified in the sample.\n")

cat("There are", sum(data$Family_ID == "indo1319", na.rm = TRUE), "individuals who list Tok Pisin as their language.\n")

cat("There are", sum(data$assigned_language_code == "tpi", na.rm = TRUE), "individuals who list Tok Pisin as their language.\n")
cat("There are", sum(data$Family_ID == "aust1307", na.rm = TRUE), "individuals who list an Austronesian language.\n")

write.csv(data, 'processed_data/genetic_samples.csv', row.names = FALSE, na = "")

