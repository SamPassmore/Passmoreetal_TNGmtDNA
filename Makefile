## This Makefile runs the analysis for Passmore et al. (2025)

## Parameters & File Paths

# How many trees to analyse
N_TREE=1000
# How many permutations for our permutation test
PERMUT=5

# Where the phylogenetic results are saved
BASE_DIR=/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/
POSTERIOR=$(BASE_DIR)/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.trees
TREE=$(BASE_DIR)/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips.trees

# Location of processed metadata
METADATA=processed_data/TNGData.csv

# Change for your BEAST location
## Or install here: https://www.beast2.org/
BEAST=/Applications/BEAST\ 2.7.6/bin/beast

help: ## Displays targets and descriptions
	@echo "Available make targets:"
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z0-9_-]+:.*?## / {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

data_processing: ## Process the genetic and geographic data for phylogenetic analysis
	Rscript processing/format_geneticdata.R
	Rscript processing/create_reference_genomes.R
	Rscript processing/distance_measures.R
	
make_phylogeny: ## Build the phylogeny using BEAST2
	# This line should run the code to build the phylogeny, although note, we did actually run the analysis on the NCI
	$(BEAST) -thread 4 analysis/20250106_LanguageSample_CBES_RelaxedClock_Constrained.xml

prepare_tree: ## Perform cleaning steps from the BEAST2 output. Some steps are manual (see commments).
	## To prepare the tree I do the following:
	### 1. Remove the burnin from the .tree file using LogCombiner.
	#### The %age of burnin is decided from Tracer, ensuring all parameters ar ESS >200
	### 2. Create the MCCT using the burnin removed file using Treecombiner
	### 3. Add the duplicate samples as contemporary polytomies in the posterior and MCCT tree using the add_duptips.R
	Rscript processing/add_duptips.R

### Run Phylogenetic Inference
# I have removed the analysis on individual languages because of memory problems
# and it is not currently discussed in the manuscript as it stands. 
individual: TNG TNG_Branch ## Run all analysis on a participants language response
Father: Father_TNG Father_Branch ## Run all analysis on a participants mother's language
Mother: Mother_TNG Mother_Branch ## Run all analysis on a participants father's language

TNG:
	echo $(TREE)
	Rscript analysis/posterior_trait_heritage_time.R \
  -t "$(TREE)" \
  -v TNG \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS

TNG_Branch:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v TNG_Branch \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS

Father_Branch:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v Branch_Family.father \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS

Father_TNG:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v TNG.father \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS

Mother_TNG:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v TNG.mother \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS

Mother_Branch:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v Branch_Family.mother \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS

All_Branches:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v Branch_Family \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS
  
TNG_Distance:
	RScript analysis/RelativeImportance.R
  
### Make Graphs
figures: ## Make figures and tables for the paper
	Rscript analysis/base_map.R
	Rscript analysis/map.R # map figure for paper
	Rscript analysis/haplogroup_table.R
	Rscript analysis/tree_plots.R
	Rscript analysis/by_languageclade.R
	# distance control
	
	