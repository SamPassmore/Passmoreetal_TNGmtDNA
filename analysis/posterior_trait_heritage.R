## Trait Heritage analysis

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
  
  # Plot functions
  library(ggplot2)
  library(ggtree)
  library(treeio)
  library(patchwork)
  
  # To parse aruguments
  library(optparse)
})

## Parse Arguments
option_list = list(
  make_option(
    c("-t", "--tree"),
    action = "store",
    default = "NCI/20240528_mtDNA_CCP/CCP-papua_m_noburnin_pruned.trees",
    type = 'character',
    help = "the location of the tree to use"
  ),
  make_option(
    c("-v", "--variable"),
    action = "store",
    default = "Branch_Family",
    type = 'character',
    help = "the variable to analyse"
  ),
  make_option(
    c("-s", "--state"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "the state within the variable to analyse"
  )
)

opt = parse_args(OptionParser(option_list=option_list))
tree_file = opt$tree
trait_name = opt$variable

# n_cores = parallel::detectCores()
n_cores = 9
n_tree = 50

## Read in the tree ####
tree = ape::read.nexus(tree_file)

tree = tree[(length(tree)-n_tree):length(tree)]

## Prepare filenames for saving
base_name = paste0("results/", trait_name, "_", tools::file_path_sans_ext(basename(tree_file)))

### Clean the tip labels to match IDs
tree = lapply(tree, function(t) {
  cleaned_tips = gsub("_.*","", t$tip.label)
  t$tip.label = cleaned_tips
  t
})

#### Read in the trait data ####
metadata = read.csv("processed_data/genetic_samples.csv")

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

## Calculate cuts based on a generation time
generation_time = 25

cat("Calculating the probability of sharing a language given the phylogeny...\n")

language_probability = pbmclapply(tree_pruned[1:n_tree], function(t)
  traitheritage::trait_heritage(t, trait = trait, generation_time = generation_time), 
  mc.preschedule = TRUE, 
  mc.cores = n_cores)

language_probability = dplyr::bind_rows(language_probability, .id = "posterior")

cat("Output saved for the probability of sharing a language given the phylogeny.\n")
write.csv(language_probability, paste0(base_name, "_languagefamily.csv"), row.names = FALSE)

# cat("Calculating the probability of sharing a language by chance (permutation tests)...\n")
# 
# permute_trait_heritage = function(tree, trait, generation_time, n_permutations = 100){
#   # tmp directory
#   dir.create("./permutations_tmp/", showWarnings = FALSE)
#   
#   # Permuted results
#   for(i in 1:n_permutations){
#     p_trait = trait
#     names(p_trait) = names(p_trait)[sample(1:length(p_trait))]
#     th = trait_heritage(tree, p_trait, generation_time)
#     write.csv(th, paste0("./permutations_tmp/", i, ".csv"))
#   }
# }
# 
# perm_prob = pbmclapply(tree_pruned[1:n_tree], function(t)
#   permute_trait_heritage(t, trait = trait, 
#                                         generation_time = generation_time, 
#                                         n_permutations = 5), 
#   mc.preschedule = TRUE, 
#   mc.cores = n_cores)
# 
# permutation_probability = lapply(perm_prob, "[[", "permuted")
# permutation_probability = dplyr::bind_rows(permutation_probability, .id = "tree_name")
# permutation_probability = permutation_probability[,-c("iteration")]

#### Distance Probability ####

metadata_dist = metadata_ss %>%
  dplyr::filter(!is.na(latitude_ofsampling) & !is.na(longitude_ofsampling))

## Calculate distance on a sphere & convert units to Kms
distance_matrix = geosphere::distm(metadata_dist[,c("longitude_ofsampling", "latitude_ofsampling")], fun = distHaversine) / 1000
dimnames(distance_matrix) = list(metadata_dist$ID, metadata_dist$ID)

tree_dist = lapply(tree, function(t) keep.tip(t, rownames(distance_matrix)))

cat("Calculating the probability of of being within 50kms of a genetic neighbour")
distance_heritage_50 = pbmclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 50
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

cat("Calculating the probability of of being within 35kms of a genetic neighbour")
distance_heritage_35 = pbmclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 35
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

cat("Calculating the probability of of being within 20kms of a genetic neighbour")
distance_heritage_20 = pbmclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 20
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_50 = dplyr::bind_rows(distance_heritage_50, .id = "posterior")
distance_heritage_35 = dplyr::bind_rows(distance_heritage_35, .id = "posterior")
distance_heritage_20 = dplyr::bind_rows(distance_heritage_20, .id = "posterior")

write.csv(permutation_probability, paste0(base_name, "_permuted_languagefamily.csv"), row.names = FALSE)
write.csv(distance_heritage_50, paste0(base_name, "_distance50.csv"), row.names = FALSE)
write.csv(distance_heritage_35, paste0(base_name, "_distance35.csv"), row.names = FALSE)
write.csv(distance_heritage_20, paste0(base_name, "_distance20.csv"), row.names = FALSE)
