## Coevolution analysis
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
    default = "./processed_data/hikingtime_h.RDS",
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
    c("-c", "--condition"),
    action = "store",
    default = 1,
    type = 'character',
    help = "the number of trees to analyse for the permutations step"
  ),
  make_option(
    c("-w", "--time"),
    action = "store",
    default = 24,
    type = 'numeric',
    help = "number of hours distance away"
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
condition = opt$condition
hours = opt$time

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

### Coevolution analysis 

cat("Calculating Coevolution....")

distance_matrix = readRDS(opt$distance)
distance_matrix = distance_matrix[rownames(distance_matrix) %in% names(trait),
                                  colnames(distance_matrix) %in% names(trait)]


tree_dist = lapply(tree_pruned[1:n_tree], function(t) keep.tip(t, rownames(distance_matrix)))

cat("Calculating the coevolution probability of of being within", hours, "hours of a genetic neighbour\n")

if(trait_name == "TNG"){
  distance_heritage_48h = mclapply(tree_dist[1:n_tree], function(t)
    traitheritage::trait_coevolution_specific(
      tree = t,
      trait = trait,
      distance_matrix = distance_matrix,
      generation_time = generation_time,
      cut_off = hours,
      condition = condition
    ), mc.preschedule = TRUE,
    mc.cores = n_cores)
} else {
  distance_heritage_48h = mclapply(tree_dist[1:n_tree], function(t)
    traitheritage::trait_coevolution(
      tree = t,
      trait = trait,
      distance_matrix = distance_matrix,
      generation_time = generation_time,
      cut_off = hours
    ), mc.preschedule = TRUE,
    mc.cores = n_cores)
}

coevolution_48 = lapply(distance_heritage_48h, "[[", "time_df")
coevolution48_conditional = data.table::rbindlist(coevolution_48, idcol = "posterior")
write.csv(coevolution48_conditional, paste0(base_name, "_coevoTNG_", hours, "h.csv"), row.names = FALSE)

rm(coevolution_48)
rm(coevolution48_conditional)
rm(distance_heritage_48h)

cat("Calculating the coevolution probability of of being within", hours, "hours of a genetic neighbour permutation\n")

if(trait_name == "TNG"){
  distance_heritage_48h_p = list()
  for(i in 1:n_tree){
    t = tree_dist[[i]]
    cat("Tree", i, "\n")
    distance_heritage_48h_p[[i]] =
      traitheritage::trait_coevolution_permutation_specific(
        tree = t,
        trait = trait,
        distance_matrix = distance_matrix,
        generation_time = generation_time,
        cut_off = hours,
        condition = condition,
        n_permutation = p_tree
      )
  }
  # distance_heritage_48h_p = mclapply(tree_dist[1:n_tree], function(t) {
  #   traitheritage::trait_coevolution_permutation_specific(
  #     tree = t,
  #     trait = trait,
  #     distance_matrix = distance_matrix,
  #     generation_time = generation_time,
  #     cut_off = hours,
  #     condition = condition,
  #     n_permutation = p_tree
  #   ) 
  # }, mc.preschedule = TRUE,
  #   mc.cores = 2)
  } else {
  # distance_heritage_48h = mclapply(tree_dist[1:n_tree], function(t)
  #   traitheritage::trait_coevolution(
  #     tree = t,
  #     trait = trait,
  #     distance_matrix = distance_matrix,
  #     generation_time = generation_time,
  #     cut_off = hours
  #   ), mc.preschedule = TRUE,
  #   mc.cores = n_cores)
}

coevolution_48 = lapply(distance_heritage_48h_p, "[[", "by_trait")
coevolution48_conditional = data.table::rbindlist(coevolution_48, idcol = "posterior")
write.csv(coevolution48_conditional, paste0(base_name, "_coevoTNG_", hours,"h_permu.csv"), row.names = FALSE)

rm(coevolution_48)
rm(coevolution48_conditional)
rm(distance_heritage_48h_p)
