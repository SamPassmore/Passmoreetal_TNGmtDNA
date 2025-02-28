library(tidyr)
library(dplyr)

dd = read.csv("/Volumes/PassmoreSSD/PapuanPast/otherdata/harley2006_table5.6.csv", skip = 1, row.names = NULL)

dd$Family[dd$Family == "East Ipili"] = "East.Ipili"

mm = as.matrix(dd[,4:ncol(dd)])
mm[lower.tri(mm)] = -99

dd[,4:ncol(dd)] = mm

TNG_Family = dd[,2:3]


tmp = dd %>% 
  pivot_longer(-one_of("Family", "TNG", "X"), names_to = "Family.cols", values_to = "signif") %>%
  filter(signif > 0) %>% 
  left_join(., TNG_Family, by = c("Family.cols" = "Family"), suffix = c("", ".cols")) %>% 
  mutate(n_comparisons = nrow(.) / 2,
         type = ifelse(TNG == TNG.cols & TNG == "Yes", "TNG", 
                       ifelse(TNG != TNG.cols, "Diff", 
                              ifelse(TNG == TNG.cols & TNG == "No", "Not", NA)))) %>% 
  group_by(type) %>% 
  summarise(n_signif = sum(signif, na.rm = T),
            n = n(),
            first(n_comparisons),
            prop = n_signif / unique(n_comparisons))
tmp
