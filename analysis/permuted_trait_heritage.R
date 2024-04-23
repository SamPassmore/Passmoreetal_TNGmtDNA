## Trait Heritage analysis

suppressPackageStartupMessages({
  # Analysis and timing packages
  library(traitheritage)
  library(ape)
  library(stringr)
  library(tictoc)
  
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
    default = "analysis/daru_n_stricprior/bayesSLC_GTR4_darun-papua_n_daru.mcct.trees",
    type = 'character',
    help = "the location of the tree to use"
  ),
  make_option(
    c("-v", "--variable"),
    action = "store",
    default = "interview_languagefamily",
    type = 'character',
    help = "the variable to analyse"
  ),
  make_option(
    c("-p", "--permutation"),
    action = "store",
    default = 100,
    type = 'numeric',
    help = "the number of permutations to complete"
  )
)

opt = parse_args(OptionParser(option_list=option_list))

tree_file = opt$tree
trait_name = opt$variable
n_permutations = opt$permutations

## Read in the tree ####
tree = read.nexus(tree_file)

### Clean the tip labels to match IDs
cleaned_tips = gsub("_.*","", tree$tip.label)
tree$tip.label = cleaned_tips

#### Read in the trait data ####
metadata = read.csv("processed_data/genetic_samples.csv")

# subset to taxa in the phylogeny
metadata_ss = metadata[metadata$ID %in% cleaned_tips,]

# how many tips are covered
sum(metadata$ID %in% cleaned_tips) / length(tree$tip.label)

# which are missing
  ## Currently the tips that are missing have been given new IDs in the 
  ## Metadata. These Ids should change when the phylogeny is re-run with new taxa labels
cleaned_tips[!cleaned_tips %in% metadata$ID]

# Isolate trait of interest
trait = metadata_ss[,trait_name]
names(trait) = metadata_ss$ID
# Missing values should be NA
trait[trait == ""] = NA
trait = trait[!is.na(trait)]

# Remove Indo-European language
trait = trait[trait != "Indo-European"]

table(trait)

## Prune the tree to match the tips
tree_pruned = keep.tip(tree, names(trait))

## Calculate cuts based on a generation time
generation_time = 25
n_cuts = floor(max(ape::node.depth.edgelength(tree_pruned)) / generation_time)

# Time how long it takes to calculate the true probability and the permuted probability 
tic()
th = traitheritage::permute_trait_heritage(tree = tree_pruned,
                                      trait = trait,
                                      n_generations = n_cuts)
toc()

## Collect extra data for plotting
plot_df = data.frame(label = names(trait),
                     trait = trait)
tree_data = full_join(tree_pruned, plot_df, by = "label")

saveRDS(list(th, tree_data), file = paste0("results/", trait_name, "_", tools::file_path_sans_ext(basename(tree_file)), ".RDS"))
