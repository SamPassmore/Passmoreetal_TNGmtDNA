library(readxl)
library(stringr)
library(dplyr)

glottolog = read.csv('https://raw.githubusercontent.com/glottolog/glottolog-cldf/master/cldf/languages.csv')

unpublished = readxl::read_xlsx("data/Languages_Nbr_1508.xlsx",
                         sheet = 1, 
                         skip = 1) 
unpublished = unpublished[,3:4]  
unpublished = unpublished[complete.cases(unpublished),]
unpublished$published = "unpublished"

published = readxl::read_xlsx("data/Languages_Nbr_1508.xlsx",
                              sheet = 2, 
                              skip = 0)
published$published = "published"

data = rbind(unpublished, published)
language_codes = str_extract_all(data$Language, "(?<=\\[).+?(?=\\])") %>% 
  plyr::ldply(., rbind)  # make the result into a matrix


data = cbind(data, language_codes)
colnames(data)[4:5] = c("Language.1", "Language.2")

## Join glottolog
data = left_join(data, glottolog, c("Language.1" = "ISO639P3code"))

# Language statistics
family_summary = table(data$Family_ID)
family_df = data.frame(family_summary)
family_df = left_join(family_df, glottolog[,c("ID", "Name")], by = c("Var1" = "ID"))

write.csv(family_df, file = "processed_data/languagefamily_count.csv", row.names = FALSE)
write.csv(data, file = "processed_data/glottolog_languagedata.csv", row.names = FALSE)
