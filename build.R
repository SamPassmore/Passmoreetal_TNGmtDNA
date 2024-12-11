## This script will build the analyses

## Install BiocManager & used packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("treeio")

# 1. Format the data received in its raw format
system("Rscript processing/format_geneticdata.R")
system("Rscript processing/create_reference_genomes.R")

# 2. Pair to glottolog
system("Rscript processing/pair_toglottolog2.R")

#. Build phylogeny

# 4. Calculate TMRCA
## For M Haplotypes
### For Men 
### For Women
## For N Haplotypes
### For Men
### For Women
