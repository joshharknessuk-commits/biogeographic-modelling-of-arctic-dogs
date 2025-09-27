#######################################################
# DEC and DEC+J
#######################################################
library(BioGeoBEARS)

# Define has.singles() if missing
if (!exists("has.singles")) {
  has.singles <- function(phy) {
    pc <- rowsum(rep(1L, nrow(phy$edge)), phy$edge[,1])
    any(pc == 1L)
  }
}


# ---- Inputs ----
trfn = "FcC_supermatrix.coyote1.rooted.fixed.tre"
geogfn = "biogeobears_input.data"
max_range_size = 2

#######################################################
# Run DEC
#######################################################
dec_run <- define_BioGeoBEARS_run()
dec_run$trfn = trfn
dec_run$geogfn = geogfn
dec_run$max_range_size = max_range_size
dec_run$min_branchlength = 0.000001
dec_run$include_null_range = TRUE
dec_run$on_NaN_error = -1e50
dec_run$speedup = TRUE
dec_run$use_optimx = "GenSA"
dec_run$num_cores_to_use = 1
dec_run$force_sparse = FALSE

dec_run <- readfiles_BioGeoBEARS_run(dec_run)
dec_run$return_condlikes_table = TRUE
dec_run$calc_TTL_loglike_from_condlikes_table = TRUE
dec_run$calc_ancprobs = TRUE

check_BioGeoBEARS_run(dec_run)

resDEC <- bears_optim_run(dec_run)
save(resDEC, file="DEC_results.Rdata")

#######################################################
# Run DEC+J
#######################################################
decj_run <- define_BioGeoBEARS_run()
decj_run$trfn = trfn
decj_run$geogfn = geogfn
decj_run$max_range_size = max_range_size
decj_run$min_branchlength = 0.000001
decj_run$include_null_range = TRUE
decj_run$on_NaN_error = -1e50
decj_run$speedup = TRUE
decj_run$use_optimx = "GenSA"
decj_run$num_cores_to_use = 1
decj_run$force_sparse = FALSE

decj_run <- readfiles_BioGeoBEARS_run(decj_run)
decj_run$return_condlikes_table = TRUE
decj_run$calc_TTL_loglike_from_condlikes_table = TRUE
decj_run$calc_ancprobs = TRUE

# Start values from DEC
decj_run$BioGeoBEARS_model_object@params_table["d","init"] = resDEC$outputs@params_table["d","est"]
decj_run$BioGeoBEARS_model_object@params_table["d","est"]  = resDEC$outputs@params_table["d","est"]
decj_run$BioGeoBEARS_model_object@params_table["e","init"] = resDEC$outputs@params_table["e","est"]
decj_run$BioGeoBEARS_model_object@params_table["e","est"]  = resDEC$outputs@params_table["e","est"]

# Add founder-event speciation
decj_run$BioGeoBEARS_model_object@params_table["j","type"] = "free"
decj_run$BioGeoBEARS_model_object@params_table["j","init"] = 0.0001
decj_run$BioGeoBEARS_model_object@params_table["j","est"]  = 0.0001
decj_run$BioGeoBEARS_model_object@params_table["j","min"]  = 1e-5
decj_run$BioGeoBEARS_model_object@params_table["j","max"]  = 0.99999

check_BioGeoBEARS_run(decj_run)

resDECj <- bears_optim_run(decj_run)
save(resDECj, file="DECJ_results.Rdata")

