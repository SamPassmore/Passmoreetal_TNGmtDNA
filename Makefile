
N_TREE=250
PERMUT=5
BASE_DIR=/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/
POSTERIOR=$(BASE_DIR)/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.trees
TREE=/Volumes/PassmoreSSD/PapuanPast/20250106_LanguageSample_CBES_Relaxed_Constrained/20250106_LanguageSample_CBES_RelaxedClock_Constrained-all_data.noburnin.duptips.trees
METADATA=processed_data/TNGData.csv

data_processing:
	Rscript processing/format_geneticdata.R
	Rscript processing/create_reference_genomes.R


prepare_tree: 
	## To prepare the tree I do the following:
	### 1. Remove the burnin from the .tree file using LogCombiner.
	#### The %age of burnin is decided from Tracer, ensuring all parameters ar ESS >200
	### 2. Create the MCCT using the burnin removed file using Treecombiner
	### 3. Add the duplicate samples as contemporary polytomies in the posterior and MCCT tree using the add_duptips.R
	Rscript processing/add_duptips.R

individual: TNG TNG_Branch interview_language
Father: Father_TNG Father_Branch father_language
Mother: Mother_TNG Mother_Branch mother_language

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

interview_language:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v interview_language \
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

father_language:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v father_language \
  -n 100 \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS

mother_language:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v mother_language \
  -n 100 \
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