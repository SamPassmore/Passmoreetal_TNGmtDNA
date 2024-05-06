## This script runs the stable clades analysis

suppressPackageStartupMessages({
  # Analysis and timing packages
  library(traitheritage)
  library(ape)
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
    default = "Branch_Family",
    type = 'character',
    help = "the variable to analyse"
  ),
  make_option(
    c("-n", "--ntree"),
    action = "store",
    default = 300,
    type = 'numeric',
    help = "the number of posterior trees to analyse"
  )
)

opt = parse_args(OptionParser(option_list=option_list))
tree_file = opt$tree
trait_name = opt$variable
n_tree = opt$ntree

#### Read in Data ####

tree = ape::read.nexus(tree_file)
# reduce to desired posterior length 
# I assume we want the tail sample of the MCMC chain
tree = tree[(length(tree)-n_tree):length(tree)]

## Metadata 
metadata = read.csv("processed_data/genetic_samples.csv")

# clean IDs
metadata$ID = gsub("-.*","", metadata$ID)

# subset to taxa in the phylogeny
metadata_ss = metadata[metadata$ID %in% tree[[1]]$tip.label,]

# extract the trait of interest
trait = metadata_ss[,trait_name]
trait[trait == ""] = NA
names(trait) = metadata_ss$ID

# remove missing traits
trait = trait[!is.na(trait)]

# remove austronesian
trait = trait[trait != "Austronesian"]

# match traits and tips
tree = keep.tip.multiPhylo(tree, names(trait))

# how many tips are covered
covered = sum(metadata$ID %in% tree[[1]]$tip.label) / length(tree[[1]]$tip.label)
cat("Matched", round(covered, 2) * 100, "% of the tips to the metadata.")

## Identify the time since MRCA clades for the posterior
fit = tmrca(tree, trait)

## Identify the most stable clades
stable_clades = fit %>% 
  group_by(tips) %>% 
  summarise(clade_count = length(TREE),
            clade_proportion = clade_count / length(tree),
            mean_Date = mean(Date),
            mode_node = mode(factor(Node)),
            state = first(state)) %>%
  filter(clade_proportion > 0.5)

stable_clades$n_tips = sapply(stable_clades$tips, function(x) length(unlist(strsplit(x, "; "))))

#### Saving the Output ####
tree_type = unlist(strsplit(basename(tree_file), "-"))[1]

f1 = paste0("results/", tree_type, "_stableclades_summary.csv")
f2 = paste0("results/", tree_type, "_stableclades_all.csv")

write.csv(stable_clades, f1)
write.csv(fit, f2)


