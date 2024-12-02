
N_TREE=500
PERMUT=3

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