
N_TREE=100
PERMUT=5
TREE=/Volumes/PassmoreSSD/PapuanPast/20241216_CBESbTest_LanguageSample/20241112_CBESbTest_LanguageSample-all_data.noburnin.duptips.trees
METADATA=processed_data/TNGData.csv

data_processing:
	Rscript processing/format_geneticdata.R
	Rscript processing/create_reference_genomes.R


prepare_tree: 
	Rscript analysis/add_duptips.R

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
  -v TNG_Branch.father \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS
  
Mother_Branch:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t $(TREE) \
  -v TNG_Branch.mother \
  -n $(N_TREE) \
  -m "$(METADATA)" \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS