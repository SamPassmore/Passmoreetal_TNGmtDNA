suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(pbmcapply)
  
  # To parse aruguments
  library(optparse)
})

## Parse Arguments
option_list = list(
  make_option(
    c("-t", "--tree"),
    action = "store",
    default = "NCI/20240618_R_N/CBES_bModelTest-papua_n.pruned.noburnin.trees",
    type = 'character',
    help = "the location of the tree to use"
  ),
  make_option(
    c("-n", "--ntree"),
    action = "store",
    default = 50,
    type = 'numeric',
    help = "the numer of trees to analyse"
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
    default = "language_metadata.csv",
    type = 'character',
    help = "the metadata file"
  ),
  make_option(
    c("-s", "--seed"),
    action = "store",
    default = 335677,
    type = 'numeric',
    help = "the metadata file"
  )
)

opt = parse_args(OptionParser(option_list=option_list))

set.seed(opt$seed)

tree_file = opt$tree
trait_name = opt$variable

n_cores = parallel::detectCores()
n_tree = opt$ntree

if (!dir.exists("results/")) {dir.create("results/")}
base_name = paste0("results/", trait_name, "_", tools::file_path_sans_ext(basename(tree_file)))


## Read in the tree ####
tree_posterior = ape::read.nexus(tree_file)
# Sub-sample posterior
tree = tree[sample(seq_along(tree_posterior), n_tree)]

# subset to the posterior sample
m_file = opt$metadata
metadata = read.csv(m_file)

metadata_ss = metadata[metadata$ID %in% tree_posterior[[1]]$tip.label, ]
trait = metadata_ss[,trait_name]
names(trait) = metadata_ss$ID

w_b = pbmcapply::pbmclapply(
  tree_posterior, function(t){
    traitheritage::within_between(tree = t, generation_time = 25, trait = trait)    
  }
)

saveRDS(w_b, paste0(base_name, "_withinbetween.RDS"))
