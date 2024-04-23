# Make XML files

library(babette)

# Data
fasta_filename = "code/processed_data/papua_m.fasta"

rates = c(2.62E-08, 2.77E-08)
freq_equilibrium = c("estimated", "empirical")
clock = c("strict", "relaxed")

settings = expand.grid(rates = rates,
                       clock = clock,
                       freq_equilibrium = freq_equilibrium
                       )
                       

## Always use the same MCMC
mcmc = beautier::create_mcmc(chain_length = 50000000 / 1000, 
                             store_every = 1000,
                             pre_burnin = 12500000 / 1000
                             )


beast2_options <- create_beast2_options()

inference_models = list()
out = list()
# for(i in 1:nrow(settings)){
for(i in 1){
  
  # Make the molecular clock model
  if(settings$clock[i] == "strict"){
    clock = create_strict_clock_model(
      clock_rate_param = settings$rates[i]
    )
  } 
  if(settings$clock[i] == "relaxed"){
    clock = create_clock_model_rln(mean_clock_rate = settings$rates[i])
  }
  
  # Make the site model
  site = create_site_model_hky(freq_equilibrium = settings$freq_equilibrium[i])
  
  # Make the tree model
  tree = create_tree_prior_cbs() # Coalescent Bayesian Skyline
  
  # Make the model 
  inference_models[[i]] = create_inference_model(
    site_model = site,
    clock_model = clock,
    tree = tree,
    mcmc = mcmc,
    beauti_options = beautier::create_beauti_options()
  )
  
  
  
  out[[i]] = bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_models[[i]],
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}

prepare_file_creation(inference_models[[1]], beast2_options = create_beast2_options())
