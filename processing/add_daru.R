suppressPackageStartupMessages({
  library(ape)
  library(stringr)
  library(Biostrings)
  library(readxl)
  library(dplyr)
})

snake_case = function(nms) tolower(gsub("[. ]+", "_", nms))
dmd = read_xlsx("data/Daru_Participant_Information_LanguageGroupings.xlsx", sheet = 1, .name_repair = snake_case)

## Read in Fasta data
m_haplotypes = read.FASTA("processed_data/papua_m.fasta")
n_haplotypes = read.FASTA("processed_data/papua_n.fasta")

daru_data = read.FASTA("data/haplogrep3_MSA.fasta_daru.txt")

## Identify Daru Duplicates
genetic_data = read.csv('processed_data/genetic_samples.csv')

n_haplotypes

sort(table(str_extract(string = genetic_data$ID, pattern = "5[0-9]{2}")))

genetic_data$ID[str_detect(genetic_data$ID, pattern = "547")]

genetic_data$ID[c(885,944)]

m_haplotypes[c("HGDP00547", "OV0547")]
dist(c(n_haplotypes[c("C001KN2-PNG540")], daru_data["540_S2"]))

## Identify haplotypes for daru data
ids = as.numeric(sapply(str_split(names(daru_data), "_"), "[[", 1))
haplogroup = sapply(str_split(names(daru_data), "_"), "[[", 2)
m_haplogroup_labels = c("M", "G", "D", "Q", "C", "Z", "E")
n_haplogroup_labels = c("N", "A", "X", "R", "B", "U", "O", "S", "I", "W", "Y", "P", "HV", "H", "V", "J", "T", "F", "K")
haplogroup_MN = ifelse(substr(haplogroup, 1, 1) %in% m_haplogroup_labels, "M",
                                    ifelse(substr(haplogroup, 1, 1) %in% n_haplogroup_labels, "N",
                                           NA))
df_join = data.frame(id = ids, haplogroup, haplogroup_MN)

dmd = left_join(dmd, df_join, by = c("sampling_id" = "id"))

dmd$sampling_id = paste0(dmd$sampling_id, "_", dmd$haplogroup)

daru_m = daru_data[haplogroup_MN == "M"]
daru_n = daru_data[haplogroup_MN == "N"]

## Add Daru M & N to main samples
m_haplotypes_wDaru = c(m_haplotypes, daru_m)
n_haplotypes_wDaru = c(n_haplotypes, daru_n)

names(m_haplotypes_wDaru[duplicated(m_haplotypes_wDaru) | duplicated(m_haplotypes_wDaru, fromLast = T)])

## Remove duplicates from Papuan past & Daru samples
### I want to keep ANU samples and Papuan Past metadata 
m_dups = m_haplotypes_wDaru[duplicated(m_haplotypes_wDaru)]

for(i in 1:length(m_dups)){
  if(isTRUE(all.equal.DNAbin(m_dups["577_Q1a1"], m_dups[i])))
    print(names(m_dups)[i])
}

dd = dist.dna(m_dups)
dist_zero = which(as.matrix(dd) == 0, arr.ind = TRUE)
dist_zero[dist_zero[,1] != dist_zero[,2],]

dist.dna(m_dups[c("577_Q1a1", "577_Q1a1")], model = "N")
dist.dna(m_dups[c("577_Q1a1", 'HRE0008')], model = "N")

View(as.matrix(dd))

all.equal.DNAbin(m_dups["519_Q1"], m_dups["531_Q3a1"])
dist.dna(m_dups[c("519_Q1", "531_Q3a1")], model = "N")

write.dna(x = m_haplotypes_wDaru, file = 'processed_data/papua_m_daru.fasta', format = "fasta")
write.dna(x = n_haplotypes_wDaru, file = 'processed_data/papua_n_daru.fasta', format = "fasta")

## Save haplogroup data to metadata 
write.csv(dmd, file = "processed_data/daru_metadata.csv", row.names = FALSE)
