# ---- Run DIVALIKE+J fresh (8 cores, no strat) ----
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(BioGeoBEARS)
})

# If your environment still complains about has.singles(), keep this small shim:
if (!exists("has.singles")) {
  has.singles <- function(phy) {
    pc <- rowsum(rep(1L, nrow(phy$edge)), phy$edge[,1])
    any(pc == 1L)
  }
}

# ---- Inputs ----
treefile <- "FcC_supermatrix.coyote1.rooted.fixed.tre"  # your rooted tree
geogfile <- "biogeobears_input_2.geog"                            # single-letter geography file
out_dir  <- getwd()
max_range_size <- 2
cores <- 8

# ---- Build a fresh run object ----
runobj <- define_BioGeoBEARS_run()
runobj$trfn <- treefile
runobj$geogfn <- geogfile
runobj$max_range_size <- max_range_size
runobj$min_branchlength <- 1e-6
runobj$include_null_range <- TRUE

runobj$on_NaN_error <- -1e50
runobj$speedup <- TRUE
runobj$use_optimx <- TRUE
runobj$num_cores_to_use <- cores
runobj$force_sparse <- FALSE

# Load inputs & defaults
runobj <- readfiles_BioGeoBEARS_run(runobj)

# Likelihood / ancestral states
runobj$return_condlikes_table <- TRUE
runobj$calc_TTL_loglike_from_condlikes_table <- TRUE
runobj$calc_ancprobs <- TRUE

# ---- DIVALIKE+J model settings ----
M <- runobj$BioGeoBEARS_model_object@params_table

# No subset sympatry; widespread vicariance allowed
M["s","type"]  <- "fixed"; M["s","init"] <- 0; M["s","est"] <- 0
M["ysv","type"]<- "2-j"
M["ys","type"] <- "ysv*1/2"
M["y","type"]  <- "ysv*1/2"
M["v","type"]  <- "ysv*1/2"

# Widespread vicariance weight
M["mx01v","type"] <- "fixed"; M["mx01v","init"] <- 0.5; M["mx01v","est"] <- 0.5

# Enable founder-event speciation (+J)
M["j","type"] <- "free"
M["j","min"]  <- 1e-5
M["j","max"]  <- 1.99999
M["j","init"] <- 0.0001
M["j","est"]  <- 0.0001

# Free d and e
M["d","init"] <- 0.01; M["d","est"] <- 0.01; M["d","type"] <- "free"
M["e","init"] <- 0.01; M["e","est"] <- 0.01; M["e","type"] <- "free"

# Attach table back to the run object
runobj$BioGeoBEARS_model_object@params_table <- M

# Apply min/max fixer
runobj <- fix_BioGeoBEARS_params_minmax(runobj)

# Sanity check
check_BioGeoBEARS_run(runobj)

# ---- Optimize model ----
resDIVALIKEJ <- bears_optim_run(runobj)

# ---- Save results ----
save(resDIVALIKEJ, file = file.path(out_dir, "DIVALIKEJ_results.Rdata"))

cat("\n--- DIVALIKE+J run complete ---\n")
cat("Saved to:", file.path(out_dir, "DIVALIKEJ_results.Rdata"), "\n")

# ---- Run DIVALIKE (no +J) on the same inputs ----
runDIV <- runobj  # copy the same run settings

# Start from a clean model table
M2 <- runDIV$BioGeoBEARS_model_object@params_table

# DIVALIKE speciation pieces (same as your +J block, but j fixed to 0)
M2["s","type"]   <- "fixed"; M2["s","init"] <- 0;    M2["s","est"] <- 0
M2["ysv","type"] <- "2-j"
M2["ys","type"]  <- "ysv*1/2"
M2["y","type"]   <- "ysv*1/2"
M2["v","type"]   <- "ysv*1/2"

# Widespread vicariance weight (keep same)
M2["mx01v","type"] <- "fixed"; M2["mx01v","init"] <- 0.5; M2["mx01v","est"] <- 0.5

# Disable founder-event speciation
M2["j","type"] <- "fixed"
M2["j","min"]  <- 0
M2["j","max"]  <- 0
M2["j","init"] <- 0
M2["j","est"]  <- 0

# Free d and e
M2["d","type"] <- "free";  M2["d","init"] <- 0.01; M2["d","est"] <- 0.01
M2["e","type"] <- "free";  M2["e","init"] <- 0.01; M2["e","est"] <- 0.01

# Attach & fix bounds
runDIV$BioGeoBEARS_model_object@params_table <- M2
runDIV <- fix_BioGeoBEARS_params_minmax(runDIV)
check_BioGeoBEARS_run(runDIV)

# Optimize
resDIVALIKE <- bears_optim_run(runDIV)

# Save
save(resDIVALIKE, file = file.path(out_dir, "DIVALIKE_results.Rdata"))
cat("\n--- DIVALIKE run complete ---\n")
cat("Saved to:", file.path(out_dir, "DIVALIKE_results.Rdata"), "\n")
