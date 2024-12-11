
N_TREE=500
PERMUT=3

data_processing:
	Rscript processing/format_geneticdata.R
	Rscript processing/create_reference_genomes.R

TNG:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t data/full_model-all.noburnin.trees \
  -v TNG \
  -n $(N_TREE) \
  -m data/TNGData.csv \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS
  
TNG_Branch:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t data/full_model-all.noburnin.trees \
  -v TNG_Branch \
  -n $(N_TREE) \
  -m data/TNGData.csv \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS
  
Father_Branch:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t data/full_model-all.noburnin.trees \
  -v TNG_Branch.father \
  -n $(N_TREE) \
  -m data/TNGData.csv \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS
  
Mother_Branch:
	Rscript analysis/posterior_trait_heritage_time.R \
  -t data/full_model-all.noburnin.trees \
  -v TNG_Branch.mother \
  -n $(N_TREE) \
  -m data/TNGData.csv \
  -p $(PERMUT) \
  -d data/hikingtime_h.RDS