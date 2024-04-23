# make hsd file

library(assertthat)
library(stringr)

# hsd Format
# 
# You can also specify your profiles in the original HaploGrep hsd format, which is a simple tab-delimited 
# file format consisting of 4 columns (ID, Range, Haplogroup and Polymorphisms).
# 
# Sample1 1-16569 H100 263G 315.1C 750G	1041G	1438G	4769G	8860G	9410G	12358G	13656C	15326G	16189C	16192T	16519C
# Sample2 1-16569 ? 73G	263G	315.1C 750G	1438G	3010A	3107C	4769G	5111T	8860G	10257T	12358G	15326G	16145A	16222T	16519C

# Genetic data processing occurs inside the Papuasphere dataset 
genetic_data = read.csv('data/genetic_samples.csv')

n_rows = nrow(genetic_data)
hsd = data.frame(ID = character(length = n_rows),
                 Haplogroup = character(length = n_rows),
                 Range = character(length = n_rows),
                 Polymorphisms = character(length = n_rows))

# ID  
hsd$ID = genetic_data$ID

# Range
## Guessing this range?
hsd$Range = "1-16569"

# Haplogroup
hsd$Haplogroup = genetic_data$haplogroup

# Polymorphisms
## Polymorphisms must be tab seperated
pm = str_remove_all(genetic_data$mutations_fromRCRS, "[,]") 
pm = str_replace_all(pm, pattern = "[ \t\n\r\v\f]", replacement = "\t")

# polymorphisms seem to have an odd structure
# remove leading letter
pm = lapply(pm, function(p){
  pp = unlist(strsplit(p, "\t"))
  
  ppp = ifelse(str_detect(pp, "[A-Z]{1}[0-9]+[A-Z]{1}"), 
               str_sub(pp, 2, nchar(pp)),
               pp)
  ppp = paste(ppp, collapse = "\t")
}) |> 
  unlist()

hsd$Polymorphisms = pm

#### Tests ####
assert_that(noNA(hsd))
assert_that(!any(str_detect(hsd$Polymorphisms, pattern = "[ \n\r\v\f]")))

write.table(
  hsd_m[1,],
  "processed_data/genetic_data_m.hsd",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  hsd_n[1,],
  "processed_data/genetic_data_n.hsd",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ape::write.nexus.data(hsd_n[,1])
