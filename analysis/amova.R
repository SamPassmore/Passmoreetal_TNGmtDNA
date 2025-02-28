library(poppr)
library(apex)
library(dplyr)
library(assertthat)
library(pegas)
library(adegenet)
library(mmod)
library(phangorn)
library(seqinr)


# Load genetic alignments
# fasta_files = list.files("processed_data/", pattern = ".fasta", full.names = TRUE)
# fasta_data = read.multiFASTA(fasta_files)
# fasta_data = read.fasta()
fasta_data = apex::read.multiFASTA("processed_data/all_data.fasta")
# And Genetic metadata
metadata = read.csv("processed_data/TNGData.csv")
## Add in duplicate metadata
TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% metadata$PhyID,] # subset to the duplicates only

metadata = rbind(metadata, TNG_dups)

# order metadata to match fasta data
metadata = metadata[match(fasta_data@labels, metadata$PhyID),]

# Check that metadata and fasta align
assert_that(sum(fasta_data@labels %in% metadata$PhyID) == 317)

# determine which language families we have data for
families = metadata %>% 
  group_by(Branch_Family) %>% 
  summarise(n = n()) %>% 
  filter(n >= 8) %>% 
  pull(Branch_Family) %>%
  na.omit()

# Reduce the genetic data to language families with sufficient coverage
keep_individuals = metadata$PhyID[metadata$Branch_Family %in% families]
fasta_data_wLanguage = fasta_data[keep_individuals]

language_strata = metadata %>%
  filter(PhyID %in% keep_individuals) %>% 
  select(interview_language, Branch_Family, TNG, haplogroup_MN)

x = assert_that(all(names(fasta_data_wLanguage) == keep_individuals))

# Make Genind object
fasta_data.gid = multidna2genind(fasta_data_wLanguage, mlst = TRUE)

# add population statification
strata(fasta_data.gid) <- language_strata
setPop(fasta_data.gid) <- ~Branch_Family


#### Pairwise population differentiation ####
diff_stats(fasta_data.gid) # this function calculates overall Nei's Gst, Hedrick's Gst and  of the dataset

nei_gst = pairwise_Gst_Nei(fasta_data.gid, linearized = FALSE) # Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- Gst)  
d_stat = pairwise_D(fasta_data.gid, linearized = FALSE, hsht_mean = "harmonic") # Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- D)  
hedrick = pairwise_Gst_Hedrick(fasta_data.gid, linearized = FALSE)

round(nei_gst, 2)
round(d_stat, 2)
round(hedrick, 2)

cor(nei_gst, d_stat)
cor(nei_gst, hedrick)
cor(d_stat, hedrick)

#### NeighbourNet ####
mtdna_dist <- dist.multidna(fasta_data_wLanguage, pool = TRUE, model = "JC69")

neighborNet(mtdna_dist)

#### AMOVA ####
mtdna_dist <- dist.multidna(fasta_data_wLanguage, pool = TRUE)


dist_mat = concatenate(fasta_data_wLanguage)
mtdna_dist <- dist.dna(dist_mat, pairwise.deletion = TRUE, as.matrix = TRUE, model = "raw")
mtdna_dist <- dist.gene(dist_mat)
as.matrix(mtdna_dist)[1:5,1:5]
# All M and N haplotypes have distances of NA
all(is.nan(as.matrix(mtdna_dist)[language_strata$haplogroup_MN == "M", language_strata$haplogroup_MN == "N"]))
any(is.nan(as.matrix(mtdna_dist)[language_strata$haplogroup_MN == "M", language_strata$haplogroup_MN == "M"]))
any(is.nan(as.matrix(mtdna_dist)[language_strata$haplogroup_MN == "N", language_strata$haplogroup_MN == "N"]))
any(is.nan(as.matrix(mtdna_dist)))

amova_output = amova(mtdna_dist ~  TNG / Branch_Family / interview_language, data = strata(fasta_data.gid), nperm = 100)
poppr.amova(fasta_data.gid, ~TNG/Branch_Family/interview_language, clonecorrect = TRUE)

