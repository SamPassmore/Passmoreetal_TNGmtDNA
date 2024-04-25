## This script runs the stable clades analysis

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
    default = "NCI/mtDNA_M/mtDNA_M-papua_m.noburnin.trees",
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

tree = ape::read.nexus(tree_file)

tree = tree[(length(tree)-n_tree):length(tree)]
