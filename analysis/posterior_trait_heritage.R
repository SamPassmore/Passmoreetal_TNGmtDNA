## Trait Heritage analysis

suppressPackageStartupMessages({
  # Analysis and timing packages
  library(traitheritage)
  library(ape)
  library(stringr)
  library(tictoc)
  library(parallel)
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
    default = "NCI/mtDNA_N/mtDNA_N-papua_n_daru_pruned.trees",
    type = 'character',
    help = "the location of the tree to use"
  ),
  make_option(
    c("-v", "--variable"),
    action = "store",
    default = "interview_languagefamily",
    type = 'character',
    help = "the variable to analyse"
  )
)

opt = parse_args(OptionParser(option_list=option_list))
tree_file = opt$tree
trait_name = opt$variable

# n_cores = parallel::detectCores()
n_cores = 8
n_tree = 100

## Read in the tree ####
tree = ape::read.tree(tree_file)

tree = tree[(length(tree)-n_tree):length(tree)]

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
cat("Matched", round(covered, 2) * 100, "% of the tips to the metadata.")

# which are missing?

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

# Time how long it takes to calculate the true probability and the permuted probability 
language_probability = parallel::mclapply(tree_pruned[1:n_tree], function(t)
  traitheritage::trait_heritage(t, trait = trait, generation_time = generation_time), 
  mc.preschedule = TRUE, 
  mc.cores = n_cores)

language_probability = dplyr::bind_rows(language_probability, .id = "posterior")

perm_prob = parallel::mclapply(tree_pruned[1:n_tree], function(t)
  traitheritage::permute_trait_heritage(t, trait = trait, generation_time = generation_time, n_permutations = 10), 
  mc.preschedule = TRUE, 
  mc.cores = n_cores)

permutation_probability = lapply(perm_prob, "[[", "permuted")
permutation_probability = dplyr::bind_rows(permutation_probability, .id = "tree_name")
permutation_probability = permutation_probability[,-c("iteration")]

#### Distance Probability ####
metadata_dist = metadata_ss %>% 
  dplyr::filter(!is.na(latitude_ofsampling) & !is.na(longitude_ofsampling))

## Calculate distance on a sphere & convert units to Kms
distance_matrix = geosphere::distm(metadata_dist[,c("longitude_ofsampling", "latitude_ofsampling")], fun = distHaversine) / 1000
dimnames(distance_matrix) = list(metadata_dist$ID, metadata_dist$ID)

tree_dist = lapply(tree, function(t) keep.tip(t, rownames(distance_matrix)))

distance_heritage_50 = parallel::mclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 50
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_35 = parallel::mclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 35
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_20 = parallel::mclapply(tree_dist[1:n_tree], function(t)
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

base_name = paste0("results/", trait_name, "_", tools::file_path_sans_ext(basename(tree_file)))

write.csv(language_probability, paste0(base_name, "_languagefamily.csv"), row.names = FALSE)
write.csv(permutation_probability, paste0(base_name, "_permuted_languagefamily.csv"), row.names = FALSE)
write.csv(distance_heritage_50, paste0(base_name, "_distance50.csv"), row.names = FALSE)
write.csv(distance_heritage_35, paste0(base_name, "_distance35.csv"), row.names = FALSE)
write.csv(distance_heritage_20, paste0(base_name, "_distance20.csv"), row.names = FALSE)
