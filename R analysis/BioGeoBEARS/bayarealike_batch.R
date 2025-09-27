#######################################################
#######################################################
# BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
#######################################################
#######################################################

library(ape)
library(phangorn)
library(BioGeoBEARS)

# ---- Inputs ----
trfn = "FcC_supermatrix.coyote1.rooted.fixed.tre"
geogfn = "biogeobears_input_2.geog"
max_range_size = 2

#######################################################
# Run BAYAREALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.0

# Linkage adjustments
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Exact copying for sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# --- Automatically fix any out-of-bounds inits ---
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "BAYAREALIKE_results.Rdata"
if (runslow) {
    resBAYAREALIKE = bears_optim_run(BioGeoBEARS_run_object)
    save(resBAYAREALIKE, file=resfn)
} else {
    load(resfn)
}

#######################################################
# Run BAYAREALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1e50
BioGeoBEARS_run_object$speedup = TRUE
BioGeoBEARS_run_object$use_optimx = "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

# Start values from BAYAREALIKE
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Allow jump dispersal
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001

# Linkage adjustments
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Exact copying for sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Safe min/max
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

# --- Automatically fix any out-of-bounds inits ---
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "BAYAREALIKEJ_results.Rdata"
if (runslow) {
    resBAYAREALIKEj = bears_optim_run(BioGeoBEARS_run_object)
    save(resBAYAREALIKEj, file=resfn)
} else {
    load(resfn)
}

