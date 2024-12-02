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
  library(data.table)
  
  # To parse aruguments
  library(optparse)
})

## Parse Arguments
option_list = list(
  make_option(
    c("-t", "--tree"),
    action = "store",
    default = "./data/full_model-all.noburnin.trees",
    type = 'character',
    help = "the location of the tree to use"
  ),
  make_option(
    c("-n", "--ntree"),
    action = "store",
    default = 3,
    type = 'numeric',
    help = "the numer of trees to analyse"
  ),
  make_option(
    c("-v", "--variable"),
    action = "store",
    default = "TNG_Branch",
    type = 'character',
    help = "the variable to analyse"
  ),
  make_option(
    c("-m", "--metadata"),
    action = "store",
    default = "data/TNGBranchData.csv",
    type = 'character',
    help = "the metadata file"
  ),
  make_option(
    c("-d", "--distance"),
    action = "store",
    default = "./data/hikingtime_h.RDS",
    type = 'character',
    help = "the distance between taxa file"
  ),
  make_option(
    c("-p", "--ptree"),
    action = "store",
    default = 3,
    type = 'numeric',
    help = "the number of trees to analyse for the permutations step"
  ),
  make_option(
    c("-s", "--seed"),
    action = "store",
    default = 93747,
    type = 'numeric',
    help = "the random seed"
  )
)

opt = parse_args(OptionParser(option_list=option_list))
tree_file = opt$tree
trait_name = opt$variable

set.seed(opt$seed)

#n_cores = parallel::detectCores()
n_cores = 4
n_tree = opt$ntree

## Read in the tree ####
trees = ape::read.nexus(tree_file)
# sample posterior
idx = sample(1:length(trees), size = n_tree)
trees = trees[idx]

## Prepare filenames for saving
if (!dir.exists("results/")) {dir.create("results/")}
base_name = paste0("results/", trait_name, "_", tools::file_path_sans_ext(basename(tree_file)))

# ### Clean the tip labels to match IDs
trees = lapply(trees, function(t) {
  cleaned_tips = gsub("_.*","", t$tip.label)
  t$tip.label = cleaned_tips
  t
})

#### Read in the trait data ####
m_file = opt$metadata
metadata = read.csv(m_file)

# clean IDs
# metadata$ID = gsub("-.*","", metadata$ID)

# subset to taxa in the phylogeny
metadata_ss = metadata[metadata$ID %in% trees[[1]]$tip.label,]

# how many tips are covered
covered = sum(metadata$ID %in% trees[[1]]$tip.label) / nrow(metadata)
cat("Matched", round(covered, 2) * 100, "% of the tips to the metadata.\n")

# Isolate trait of interest
trait = metadata_ss[,trait_name]
names(trait) = metadata_ss$ID
# Missing values should be NA
trait[trait == ""] = NA
trait = trait[!is.na(trait)]

## Prune the tree to match the tips
tree_pruned = lapply(trees, function(t) keep.tip(t, names(trait)))
class(tree_pruned) = "multiPhylo"

## Calculate cuts based on a generation time
generation_time = 25

cat("Calculating the probability of sharing a language given the phylogeny...\n")

language_probability = mclapply(tree_pruned[1:n_tree], function(t)
  traitheritage::trait_heritage(t, trait = trait, generation_time = generation_time),
  mc.preschedule = TRUE,
  mc.cores = n_cores)

languageprobability_bytrait = lapply(language_probability, function(x) x$by_trait) %>%
  data.table::rbindlist(., idcol = "posterior")

languageprobability_summary = lapply(language_probability, function(x) x$summary) %>%
  data.table::rbindlist(., idcol = "posterior")

cat("Output saved for the probability of sharing a language given the phylogeny.\n")
write.csv(languageprobability_bytrait, paste0(base_name, "_probability_bytrait.csv"), row.names = FALSE)
write.csv(languageprobability_summary, paste0(base_name, "_probability_summary.csv"), row.names = FALSE)

rm(language_probability)

cat("Calculating the probability of sharing a language by chance (permutation tests)...\n")

perm_prob = mclapply(tree_pruned[1:n_tree], function(t)
  permute_trait_heritage(t, trait = trait,
                                        generation_time = generation_time,
                                        n_permutations = 5),
  mc.preschedule = TRUE,
  mc.cores = n_cores)

permutationprobability_bytrait = lapply(perm_prob, "[[", "by_trait")
permutationprobability_summary = lapply(perm_prob, "[[", "summary")

permutationprobability_bytrait = data.table::rbindlist(permutationprobability_bytrait, idcol = "posterior")
permutationprobability_summary = data.table::rbindlist(permutationprobability_summary, idcol = "posterior")

write.csv(permutationprobability_bytrait, paste0(base_name, "_permutedbytrait.csv"), row.names = FALSE)
write.csv(permutationprobability_summary, paste0(base_name, "_permutedsummary.csv"), row.names = FALSE)

rm(perm_prob)

#### Distance Probability ####

distance_matrix = readRDS(opt$distance)
distance_matrix = distance_matrix[rownames(distance_matrix) %in% names(trait),
                                  colnames(distance_matrix) %in% names(trait)]


tree_dist = lapply(tree_pruned[1:n_tree], function(t) keep.tip(t, rownames(distance_matrix)))


cat("Calculating the probability of of being within 48 hours of a genetic neighbour\n")
distance_heritage_48h = mclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 48
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_48h_bytrait = lapply(distance_heritage_48h, "[[", "by_trait")
distance_heritage_48h_summary = lapply(distance_heritage_48h, "[[", "summary")

distance_heritage_48h_bytrait = data.table::rbindlist(distance_heritage_48h_bytrait, idcol = "posterior")
distance_heritage_48h_summary = data.table::rbindlist(distance_heritage_48h_summary, idcol = "posterior")

write.csv(distance_heritage_48h_bytrait, paste0(base_name, "_48hbytrait.csv"), row.names = FALSE)
write.csv(distance_heritage_48h_bytrait, paste0(base_name, "_48hsummary.csv"), row.names = FALSE)

rm(distance_heritage_48h)

cat("Calculating the probability of of being within 24hours of a genetic neighbour\n")
distance_heritage_24h = mclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 24
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_24h_bytrait = lapply(distance_heritage_24h, "[[", "by_trait")
distance_heritage_24h_summary = lapply(distance_heritage_24h, "[[", "summary")

distance_heritage_24h_bytrait = data.table::rbindlist(distance_heritage_24h_bytrait, idcol = "posterior")
distance_heritage_24h_summary = data.table::rbindlist(distance_heritage_24h_summary, idcol = "posterior")

write.csv(distance_heritage_24h_bytrait, paste0(base_name, "_24hbytrait.csv"), row.names = FALSE)
write.csv(distance_heritage_24h_bytrait, paste0(base_name, "_24hsummary.csv"), row.names = FALSE)

rm(distance_heritage_24h)

cat("Calculating the probability of of being within 12hours walk of a genetic neighbour\n")
distance_heritage_12h = mclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 12
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_12h_bytrait = lapply(distance_heritage_12h, "[[", "by_trait")
distance_heritage_12h_summary = lapply(distance_heritage_12h, "[[", "summary")

distance_heritage_12h_bytrait = data.table::rbindlist(distance_heritage_12h_bytrait, idcol = "posterior")
distance_heritage_12h_summary = data.table::rbindlist(distance_heritage_12h_summary, idcol = "posterior")

write.csv(distance_heritage_12h_bytrait, paste0(base_name, "_12hbytrait.csv"), row.names = FALSE)
write.csv(distance_heritage_12h_bytrait, paste0(base_name, "_12hsummary.csv"), row.names = FALSE)

rm(distance_heritage_12h)

cat("Calculating the probability of of being within 6hours of a genetic neighbour\n")
distance_heritage_6h = mclapply(tree_dist[1:n_tree], function(t)
  traitheritage::distance_trait_heritage(
    tree = t,
    generation_time = generation_time,
    distance_matrix = distance_matrix,
    cut_off = 6
  ), mc.preschedule = TRUE,
  mc.cores = n_cores)

distance_heritage_6h_bytrait = lapply(distance_heritage_6h, "[[", "by_trait")
distance_heritage_6h_summary = lapply(distance_heritage_6h, "[[", "summary")

distance_heritage_6h_bytrait = data.table::rbindlist(distance_heritage_6h_bytrait, idcol = "posterior")
distance_heritage_6h_summary = data.table::rbindlist(distance_heritage_6h_summary, idcol = "posterior")

write.csv(distance_heritage_6h_bytrait, paste0(base_name, "_6hbytrait.csv"), row.names = FALSE)
write.csv(distance_heritage_6h_bytrait, paste0(base_name, "_6hsummary.csv"), row.names = FALSE)

rm(distance_heritage_6h)
