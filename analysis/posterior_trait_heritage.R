## Trait Heritage analysis
library(profvis) #uses utils::Rprof
gc()
suppressPackageStartupMessages({
  # Analysis and timing packages
  library(traitheritage)
  library(ape)
  library(stringr)
  library(tictoc)
  library(parallel)
  library(pbmcapply)
  library(geosphere)
  library(dplyr)
  
  # To parse aruguments
  library(optparse)
})

## Parse Arguments
option_list = list(
  make_option(
    c("-t", "--tree"),
    action = "store",
    default = "NCI/20240702_N_CBES_bModelTest/20240617_CBES_bModelTest-papua_n.noburnin.trees",
    type = 'character',
    help = "the location of the tree to use"
  ),
  make_option(
    c("-n", "--ntree"),
    action = "store",
    default = 10,
    type = 'numeric',
    help = "the numer of trees to analyse"
  ),
  make_option(
    c("-p", "--ptree"),
    action = "store",
    default = 3,
    type = 'numeric',
    help = "the numer of trees to analyse for the permutations step"
  ),
  make_option(
    c("-v", "--variable"),
    action = "store",
    default = "Branch_Family",
    type = 'character',
    help = "the variable to analyse"
  ),
  make_option(
    c("-m", "--metadata"),
    action = "store",
    default = "./NCI/20240618_R_N/language_metadata.csv",
    type = 'character',
    help = "the metadata file"
  ),
  make_option(
    c("-s", "--seed"),
    action = "store",
    default = 93747,
    type = 'numeric',
    help = "the metadata file"
  )
)

opt = parse_args(OptionParser(option_list=option_list))

set.seed(opt$seed)

tree_file = opt$tree
trait_name = opt$variable

# n_cores = parallel::detectCores()
n_cores = 4
n_tree = opt$ntree
p_tree = opt$ptree

## Calculate cuts based on a generation time
generation_time = 25

## Read in the tree ####
tree = ape::read.nexus(tree_file)

# take a random sample
tree = tree[sample(seq_along(tree), n_tree)]

## Prepare filenames for saving
if (!dir.exists("results/")) {dir.create("results/")}
base_name = paste0("results/", trait_name, "_", tools::file_path_sans_ext(basename(tree_file)))

### Clean the tip labels to match IDs
tree = lapply(tree, function(t) {
  cleaned_tips = gsub("_.*","", t$tip.label)
  t$tip.label = cleaned_tips
  t
})

#### Read in the trait data ####
m_file = opt$metadata
metadata = read.csv(m_file)

# clean IDs
metadata$ID = gsub("-.*","", metadata$ID)

# subset to taxa in the phylogeny
metadata_ss = metadata[metadata$ID %in% tree[[1]]$tip.label,]

# how many tips are covered
covered = sum(metadata$ID %in% tree[[1]]$tip.label) / length(tree[[1]]$tip.label)
cat("Matched", round(covered, 2) * 100, "% of the tips to the metadata.\n")

# Isolate trait of interest
trait = metadata_ss[,trait_name]
names(trait) = metadata_ss$ID
# Missing values should be NA
trait[trait == ""] = NA
trait = trait[!is.na(trait)]

## This line is specific to Language family variable. 
## Think about how to make this script independent of a particular variable
# Remove Indo-European language
trait = trait[trait != "Indo-European"]

## Prune the tree to match the tips
tree_pruned = lapply(tree, function(t) keep.tip(t, names(trait)))
class(tree_pruned) = "multiPhylo"

cat("Calculating the probability of sharing a language given the phylogeny...\n")
n_sample = sample(seq_along(tree_pruned), n_tree, replace = FALSE)
language_probability = pbmclapply(tree_pruned[1:n_tree], function(t)
  traitheritage::trait_heritage(t, trait = trait, generation_time = generation_time), 
  mc.preschedule = TRUE, 
  mc.cores = n_cores)

language_probability = dplyr::bind_rows(language_probability, .id = "posterior")

cat("Output saved for the probability of sharing a language given the phylogeny.\n")
write.csv(language_probability, paste0(base_name, "_languagefamily.csv"), row.names = FALSE)

rm(language_probability)
gc()

cat("Calculating the probability of sharing a language by chance (permutation tests)...\n")

# Subsample for the permutation step because it is so intensive
p_sample = sample(1:length(tree_pruned), p_tree, replace = FALSE)

perm_prob = pbmclapply(tree_pruned[p_sample], function(t)
  permute_trait_heritage(t, trait = trait,
                                        generation_time = generation_time,
                                        n_permutations = 5),
  mc.preschedule = TRUE,
  mc.cores = n_cores)

permutation_probability = lapply(perm_prob, "[[", "permuted")
permutation_probability = dplyr::bind_rows(permutation_probability, .id = "tree_name")
permutation_probability = permutation_probability[,-c("iteration")]

write.csv(permutation_probability, paste0(base_name, "_permuted_languagefamily.csv"), row.names = FALSE)

rm(permutation_probability)
gc()


#### Distance Probability ####

metadata_dist = metadata_ss %>%
  dplyr::filter(!is.na(latitude_hybrid) & !is.na(longitude_hybrid))

## Calculate distance on a sphere & convert units to Kms
distance_matrix = geosphere::distm(metadata_dist[,c("longitude_hybrid", "latitude_hybrid")], fun = distHaversine) / 1000
dimnames(distance_matrix) = list(metadata_dist$ID, metadata_dist$ID)

tree_dist = lapply(tree, function(t) keep.tip(t, rownames(distance_matrix)))

cat("Calculating the probability of of being within 50kms of a genetic neighbour\n")
distance_heritage_50 = pbmclapply(tree_dist[n_sample], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 50
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_50 = dplyr::bind_rows(distance_heritage_50, .id = "posterior")
write.csv(distance_heritage_50, paste0(base_name, "_distance50.csv"), row.names = FALSE)

rm(distance_heritage_50)
gc()

cat("Calculating the probability of of being within 35kms of a genetic neighbour\n")
distance_heritage_35 = pbmclapply(tree_dist[n_sample], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 35
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_35 = dplyr::bind_rows(distance_heritage_35, .id = "posterior")
write.csv(distance_heritage_35, paste0(base_name, "_distance35.csv"), row.names = FALSE)

rm(distance_heritage_35)
gc()

cat("Calculating the probability of of being within 20kms of a genetic neighbour\n")
distance_heritage_20 = pbmclapply(tree_dist[n_sample], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 20
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_20 = dplyr::bind_rows(distance_heritage_20, .id = "posterior")
write.csv(distance_heritage_20, paste0(base_name, "_distance20.csv"), row.names = FALSE)

rm(distance_heritage_20)
gc()

beepr::beep(8)
