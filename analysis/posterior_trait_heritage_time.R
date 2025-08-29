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
    default = "/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips.trees",
    type = 'character',
    help = "the location of the tree to use"
  ),
  make_option(
    c("-n", "--ntree"),
    action = "store",
    default = 2,
    type = 'numeric',
    help = "the numer of trees to analyse"
  ),
  make_option(
    c("-v", "--variable"),
    action = "store",
    default = "TNG",
    type = 'character',
    help = "the variable to analyse"
  ),
  make_option(
    c("-m", "--metadata"),
    action = "store",
    default = "./processed_data/TNGData.csv",
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
    default = 1,
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
p_tree = opt$ptree

## Read in the tree ####
trees = ape::read.nexus(tree_file)

# sample posterior
# idx = sample(1:length(trees), size = n_tree)
idx = (length(trees) - (n_tree - 1)):length(trees)
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

TNG_dups = read.csv("processed_data/genetic_duplicate_data.csv")
TNG_dups = TNG_dups[!TNG_dups$PhyID %in% metadata$PhyID,] # subset to the duplicates only

metadata = rbind(metadata, TNG_dups)

cat("There are", nrow(metadata), "taxa in the metadata")

## Remove some errant entries
# Tok Pisin langauges should be NA since it is a recently arisen langauge.
metadata = metadata[metadata$interview_language != "Tok Pisin [tpi]",]
# Baining speakers are unlikely to have come from New Guinea Mainland. This is a metadata mistake. 
metadata = metadata[metadata$Branch_Family != "Baining",]
# subset to taxa in the phylogeny
metadata_ss = metadata[metadata$ID %in% trees[[1]]$tip.label,]

# how many tips are covered
covered = sum(metadata_ss$ID %in% trees[[1]]$tip.label) / nrow(metadata_ss)
cat("Matched", round(covered, 2) * 100, "% of the tips to the metadata.\n")

# Isolate trait of interest
trait = metadata_ss[,trait_name]
names(trait) = metadata_ss$ID
# Missing values should be NA
trait[trait == ""] = NA
# NotTNG values should be NA
# trait[trait == "NotTNG"] = NA

trait = trait[!is.na(trait)]

## Prune the tree to match the tips
tree_pruned = lapply(trees, function(t) keep.tip(t, names(trait)))
class(tree_pruned) = "multiPhylo"

## Calculate cuts based on a generation time
generation_time = 25

cat("Calculating the probability of sharing a language given the phylogeny...\n")

if(trait_name == "TNG"){
  language_probability = mclapply(tree_pruned, function(t)
    traitheritage::trait_heritage_specific(t, trait = trait, generation_time = generation_time, condition = 1),
    mc.preschedule = TRUE,
    mc.cores = n_cores)
} else {
  language_probability = mclapply(tree_pruned, function(t)
    traitheritage::trait_heritage(t, trait = trait, generation_time = generation_time),
    mc.preschedule = TRUE,
    mc.cores = n_cores)
}

languageprobability_bytrait = lapply(language_probability, function(x) x$by_trait) %>%
  data.table::rbindlist(., idcol = "posterior")

write.csv(languageprobability_bytrait, paste0(base_name, "_probability_bytrait.csv"), row.names = FALSE)
rm(languageprobability_bytrait)

languageprobability_summary = lapply(language_probability, function(x) x$summary) %>%
  data.table::rbindlist(., idcol = "posterior")

write.csv(languageprobability_summary, paste0(base_name, "_probability_summary.csv"), row.names = FALSE)
rm(languageprobability_summary)

cat("Output saved for the probability of sharing a language given the phylogeny.\n")
rm(language_probability)

cat("Calculating the probability of sharing a language by chance (permutation tests)...\n")

if(trait_name == "TNG"){
  perm_prob = mclapply(tree_pruned, function(t)
    permute_trait_heritage_specific(t, trait = trait,
                                    generation_time = generation_time,
                                    n_permutations = p_tree,
                                    condition = 1),
    mc.preschedule = TRUE,
    mc.cores = n_cores)
  
} else {
  perm_prob = mclapply(tree_pruned, function(t)
    permute_trait_heritage(t, trait = trait,
                                    generation_time = generation_time,
                                    n_permutations = p_tree),
    mc.preschedule = TRUE,
    mc.cores = n_cores)
  
}

permutationprobability_bytrait = lapply(perm_prob, "[[", "by_trait")
permutationprobability_bytrait = data.table::rbindlist(permutationprobability_bytrait, idcol = "posterior")
write.csv(permutationprobability_bytrait, paste0(base_name, "_permutedbytrait.csv"), row.names = FALSE)

rm(permutationprobability_bytrait)

permutationprobability_summary = lapply(perm_prob, "[[", "summary")
permutationprobability_summary = data.table::rbindlist(permutationprobability_summary, idcol = "posterior")
write.csv(permutationprobability_summary, paste0(base_name, "_permutedsummary.csv"), row.names = FALSE)

rm(permutationprobability_summary)


rm(perm_prob)

#### Distance Probability ####

# cat("Calculating Distance....")
# 
# sort( sapply(ls(),function(x){object.size(get(x))})) 
# 
# distance_matrix = readRDS(opt$distance)
# distance_matrix = distance_matrix[rownames(distance_matrix) %in% names(trait),
#                                   colnames(distance_matrix) %in% names(trait)]
# 
# 
# tree_dist = lapply(tree_pruned[1:n_tree], function(t) keep.tip(t, rownames(distance_matrix)))
# 
# 
# cat("Calculating the probability of of being within 48 hours of a genetic neighbour\n")
# distance_heritage_48h = mclapply(tree_dist[1:n_tree], function(t)
#   traitheritage::distance_trait_heritage(
#     tree = t,
#     generation_time = generation_time,
#     distance_matrix = distance_matrix,
#     cut_off = 48
#   ), mc.preschedule = TRUE,
#   mc.cores = n_cores)
# 
# 
# distance_heritage_48h_bytrait = lapply(distance_heritage_48h, "[[", "by_trait")
# distance_heritage_48h_bytrait = data.table::rbindlist(distance_heritage_48h_bytrait, idcol = "posterior")
# write.csv(distance_heritage_48h_bytrait, paste0(base_name, "_48hbytrait.csv"), row.names = FALSE)
# 
# rm(distance_heritage_48h_bytrait)
# 
# distance_heritage_48h_summary = lapply(distance_heritage_48h, "[[", "summary")
# distance_heritage_48h_summary = data.table::rbindlist(distance_heritage_48h_summary, idcol = "posterior")
# write.csv(distance_heritage_48h_summary, paste0(base_name, "_48hsummary.csv"), row.names = FALSE)
# rm(distance_heritage_48h_summary)
# 
# sort( sapply(ls(),function(x){object.size(get(x))})) 
# 
# rm(distance_heritage_48h)
# 
# cat("Calculating the probability of of being within 24hours of a genetic neighbour\n")
# distance_heritage_24h = mclapply(tree_dist[1:n_tree], function(t)
#   traitheritage::distance_trait_heritage(
#     tree = t,
#     generation_time = generation_time,
#     distance_matrix = distance_matrix,
#     cut_off = 24
#   ), mc.preschedule = TRUE,
#   mc.cores = n_cores)
# 
# distance_heritage_24h_bytrait = lapply(distance_heritage_24h, "[[", "by_trait")
# distance_heritage_24h_bytrait = data.table::rbindlist(distance_heritage_24h_bytrait, idcol = "posterior")
# write.csv(distance_heritage_24h_bytrait, paste0(base_name, "_24hbytrait.csv"), row.names = FALSE)
# 
# rm(distance_heritage_24h_bytrait)
# 
# distance_heritage_24h_summary = lapply(distance_heritage_24h, "[[", "summary")
# distance_heritage_24h_summary = data.table::rbindlist(distance_heritage_24h_summary, idcol = "posterior")
# write.csv(distance_heritage_24h_summary, paste0(base_name, "_24hsummary.csv"), row.names = FALSE)
# 
# rm(distance_heritage_24h_summary)
# 
# sort( sapply(ls(),function(x){object.size(get(x))})) 
# 
# rm(distance_heritage_24h)
# 
# cat("Calculating the probability of of being within 12hours walk of a genetic neighbour\n")
# distance_heritage_12h = mclapply(tree_dist[1:n_tree], function(t)
#   traitheritage::distance_trait_heritage(
#     tree = t,
#     generation_time = generation_time,
#     distance_matrix = distance_matrix,
#     cut_off = 12
#   ), mc.preschedule = TRUE,
#   mc.cores = n_cores)
# 
# distance_heritage_12h_bytrait = lapply(distance_heritage_12h, "[[", "by_trait")
# distance_heritage_12h_bytrait = data.table::rbindlist(distance_heritage_12h_bytrait, idcol = "posterior")
# write.csv(distance_heritage_12h_bytrait, paste0(base_name, "_12hbytrait.csv"), row.names = FALSE)
# 
# rm(distance_heritage_12h_bytrait)
# 
# distance_heritage_12h_summary = lapply(distance_heritage_12h, "[[", "summary")
# distance_heritage_12h_summary = data.table::rbindlist(distance_heritage_12h_summary, idcol = "posterior")
# write.csv(distance_heritage_12h_summary, paste0(base_name, "_12hsummary.csv"), row.names = FALSE)
# 
# rm(distance_heritage_12h_summary)
# 
# sort( sapply(ls(),function(x){object.size(get(x))})) 
# 
# rm(distance_heritage_12h)
# 
# cat("Calculating the probability of of being within 6hours of a genetic neighbour\n")
# distance_heritage_6h = mclapply(tree_dist[1:n_tree], function(t)
#   traitheritage::distance_trait_heritage(
#     tree = t,
#     generation_time = generation_time,
#     distance_matrix = distance_matrix,
#     cut_off = 6
#   ), mc.preschedule = TRUE,
#   mc.cores = n_cores)
# 
# distance_heritage_6h_bytrait = lapply(distance_heritage_6h, "[[", "by_trait")
# distance_heritage_6h_bytrait = data.table::rbindlist(distance_heritage_6h_bytrait, idcol = "posterior")
# write.csv(distance_heritage_6h_bytrait, paste0(base_name, "_6hbytrait.csv"), row.names = FALSE)
# 
# rm(distance_heritage_6h_bytrait)
# 
# distance_heritage_6h_summary = lapply(distance_heritage_6h, "[[", "summary")
# distance_heritage_6h_summary = data.table::rbindlist(distance_heritage_6h_summary, idcol = "posterior")
# write.csv(distance_heritage_6h_summary, paste0(base_name, "_6hsummary.csv"), row.names = FALSE)
# 
# rm(distance_heritage_6h_summary)
# 
# cat("Calculating the permuted probability of of being within 6hours of a genetic neighbour\n")
# distance_heritage_6h = mclapply(tree_dist[1:n_tree], function(t)
#   traitheritage::distance_trait_heritage(
#     tree = t,
#     generation_time = generation_time,
#     distance_matrix = distance_matrix,
#     cut_off = 6
#   ), mc.preschedule = TRUE,
#   mc.cores = n_cores)
# 
# sort( sapply(ls(),function(x){object.size(get(x))})) 
# 
# rm(distance_heritage_6h)