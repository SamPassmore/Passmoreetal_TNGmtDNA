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

data_processing: ## Process the genetic data for phylogenetic analysis
	Rscript processing/format_geneticdata.R
	Rscript processing/create_reference_genomes.R

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

branch: Angan ChimbuWahgi KainantuGoroka SukiGogodala Enga Kiwaian

TNG:
	echo $(TREE)
	Rscript analysis/posterior_trait_heritage_time.R \
  -t "$(TREE)" \
  -v TNG \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -c "TNG" 
  
TNG_Coevo48:
	echo $(TREE)
	Rscript analysis/posterior_coevolution.R \
	-t "$(TREE)" \
	-v TNG \
	-n 900 \
	-m "$(METADATA)" \
	-c "1" \
	-p 3 \
	-w 48 \ 
	-d processed_data/hikingtime_h.RDS 

TNG_Coevo24:
	echo $(TREE)
	Rscript analysis/posterior_coevolution.R \
		-t "$(TREE)" \
		-v TNG \
		-n 900 \
		-m "$(METADATA)" \
		-c "1" \
		-p 3 \
		-w 24 \
		-d processed_data/hikingtime_h.RDS 
		
TNG_Coevo12:
	echo $(TREE)
	Rscript analysis/posterior_coevolution.R \
		-t "$(TREE)" \
		-v TNG \
		-n 900 \
		-m "$(METADATA)" \
		-c "1" \
		-p 3 \
		-w 12 \
		-d processed_data/hikingtime_h.RDS 

Angan:
	echo $(TREE)
	Rscript analysis/posterior_trait_heritage_time.R \
  -t "$(TREE)" \
  -v TNG_Branch \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -c "Angan" 
  
ChimbuWahgi:
	echo $(TREE)
	Rscript analysis/posterior_trait_heritage_time.R \
  -t "$(TREE)" \
  -v TNG_Branch \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -c "Chimbu-Wahgi"
  
KainantuGoroka:
	echo $(TREE)
	Rscript analysis/posterior_trait_heritage_time.R \
  -t "$(TREE)" \
  -v TNG_Branch \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -c "Kainantu-Goroka"
  
SukiGogodala:
	echo $(TREE)
	Rscript analysis/posterior_trait_heritage_time.R \
  -t "$(TREE)" \
  -v TNG_Branch \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -c "Suki-Gogodala"
  
Enga:
	echo $(TREE)
	Rscript analysis/posterior_trait_heritage_time.R \
  -t "$(TREE)" \
  -v TNG_Branch \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -c "Enga-Kewa-Huli"
  
Kiwaian:
	echo $(TREE)
	Rscript analysis/posterior_trait_heritage_time.R \
  -t "$(TREE)" \
  -v TNG_Branch \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -c "Kiwaian"

Father_TNG:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v TNG.father \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d processed_data/hikingtime_h.RDS

Mother_TNG:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v TNG.mother \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d processed_data/hikingtime_h.RDS

### Make Graphs
figures: ## Make figures and tables for the paper
	Rscript analysis/map.R # map figure
	Rscript analysis/haplogroup_table.R # make table of haplogroup frequencies
	Rscript analysis/tree_plots.R # make the tree and probability plot
	Rscript analysis/languagesubgroup_plot.R # make the subgroup family plot
	Rscript analysis/geographic_figure.R ## Make the coevolution of geography and phylogeny plot
	