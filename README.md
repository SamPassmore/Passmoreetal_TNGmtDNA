# Code for Tracing Shared Gene-Language Histories in Papua New Guinea: Signals of a Deep Antiquity in the Trans-New Guinea Languages (Passmore et al. In Prep.)
_Sam Passmore, Ray Tobler, Nicolas Brucato, Jimmy Nébni, Monika Karmin, Matthew Leavesley, Nicholas Evans, Francois-Xavier Ricaut_

This repository holds the code to replicate much of the analysis in Passmore et al. (In Prep.). The analysis procedure is defined within the Makefile and can be run using the make commands provided (see below).
The only exception is the post-processing of the mtDNA phylogeny, which was produced using LogCombiner (Drummond and Rambaut, 2007) and  and TreeAnnotator (Heled and Bouckaert, 2013). 

```
Available make targets:
  help                 Displays targets and descriptions
  data_processing      Process the genetic data for phylogenetic analysis
  make_phylogeny       Build the phylogeny using BEAST2
  prepare_tree         Perform cleaning steps from the BEAST2 output. Some steps are manual (see commments).
  individual           Run all analysis on a participants language response
  branch               Run the analysis for each sub-group of interest
  figures              Make figures and tables for the paper
```

To run this analysis, you will be required to have the following R packages installed:
```
ape
tidyr
dplyr
ggplot2
patchwork
treeio
stringr
tictoc
paralell
geosphere
data.table
optparse
forcats
kableExtra
knitr
janitor
terra
tidyterra
sf
rnaturalearth
rnaturalearthdata
ggspatial
Biostrings
pbapply
assertthat
readxl
insight
phangorn
raster
gdistance
reshape2
```
And the custom package:
`traitheritage (github.comSamPassmore/traitheritage)`

Drummond, A. J., & Rambaut, A. (2007). BEAST: Bayesian evolutionary analysis by sampling trees. BMC Evolutionary Biology, 7(1), 214. https://doi.org/10.1186/1471-2148-7-214

Heled, J., & Bouckaert, R. R. (2013). Looking for trees in the forest: Summary tree from posterior samples. BMC Evolutionary Biology, 13(1), 221. https://doi.org/10.1186/1471-2148-13-221
