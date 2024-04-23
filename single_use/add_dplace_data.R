# match d-place data
library(dplyr)
library(tidyr)

languages = read.csv('data/language_metadata.csv')

# EA010 - EA015 are variables about residence after marriage
# EA043 Major Descent pattern
# EA043 is made up from EA017, 19, and 21. 
varsof_interest = paste0("EA0", c(10:15, 17, 19, 21, 43))
xd_id = read.csv("https://raw.githubusercontent.com/D-PLACE/dplace-data/master/datasets/EA/societies.csv")


dplace = read.csv('https://raw.githubusercontent.com/D-PLACE/dplace-data/master/datasets/EA/data.csv') %>% 
  filter(var_id %in% varsof_interest) %>% 
  pivot_wider(names_from = var_id, values_from = code) %>% 
  left_join(., xd_id, by = c("soc_id" = "id"))


languages = left_join(languages, dplace, by = "glottocode")

write.csv(languages, file = "data/language_metadata.csv")

# 1: Avunulocal   8: Patrilocal   10: Virilocal
table(languages$EA010)
# 1: Wife to Husband's group 
table(languages$EA011)
# 1: Avunulocal 4: Optionally Patrilocal 8: Patrilocal
table(languages$EA012)
# 2: No pattern or Neolocal. 3: Husband to Wife's group. 9: No alternate
table(languages$EA013)
# 6: Neolocal. 9: Uxorilocal. 10: Virilocal. 11: No alternate
table(languages$EA014)
# 2: Segmented communities. 3: Agamous communities (no tendency to exogamy or endogamy) 6: Clan communities
table(languages$EA015)
# 1: None. 3: Modest Size. 4: Sibs / Clans 5: Pharatries
table(languages$EA017)
# 1: none 3: Modest size. 4: Sibs / Clans 5: Pharatries 
table(languages$EA019)
# 1: Patrilineal; 2: Duolateral. 3: Matrilineal; 7: Mixed
table(languages$EA043)
