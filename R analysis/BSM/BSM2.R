#######################################################
# Extend DEC+J BSM from 50 to 150 total replicates
#######################################################

library(BioGeoBEARS)

model_name <- "DECJ"

# ---- 1) Load DEC+J results ----
load("DECJ_results.Rdata")   # must contain resDECj
res <- resDECj               # assign to generic name used later

# ---- 2) Load your *existing* 50 BSM reps ----
old_clado <- old_ana <- NULL

# Preferred: combined full BSM output
if (file.exists("DECJ_BSM_output_full.Rdata")) {
  load("DECJ_BSM_output_full.Rdata")
  if (exists("BSM_output")) {
    old_clado <- BSM_output$RES_clado_events_tables
    old_ana   <- BSM_output$RES_ana_events_tables
  }
}

# Fallback: separate event table files
if (is.null(old_clado) && file.exists("DECJ_RES_clado_events_tables.Rdata")) {
  load("DECJ_RES_clado_events_tables.Rdata") # should create RES_clado_events_tables
  old_clado <- RES_clado_events_tables
}
if (is.null(old_ana) && file.exists("DECJ_RES_ana_events_tables.Rdata")) {
  load("DECJ_RES_ana_events_tables.Rdata")   # should create RES_ana_events_tables
  old_ana <- RES_ana_events_tables
}

# Absolute fallback: generic names
if (is.null(old_clado) && file.exists("RES_clado_events_tables.Rdata")) {
  load("RES_clado_events_tables.Rdata")
  old_clado <- RES_clado_events_tables
}
if (is.null(old_ana) && file.exists("RES_ana_events_tables.Rdata")) {
  load("RES_ana_events_tables.Rdata")
  old_ana <- RES_ana_events_tables
}

if (is.null(old_clado) || is.null(old_ana)) {
  stop("Could not find existing BSM event tables.")
}

message(sprintf("Loaded existing BSMs: %d clado maps, %d ana maps",
                length(old_clado), length(old_ana)))

# ---- 3) Prepare inputs for new BSM run ----
stochastic_mapping_inputs_list <- get_inputs_for_stochastic_mapping(res = res)

# ---- 4) Run an additional 100 BSMs ----
set.seed(54321)
BSM_output_new <- runBSM(
  res,
  stochastic_mapping_inputs_list = stochastic_mapping_inputs_list,
  maxnum_maps_to_try = 100,
  nummaps_goal = 100,
  maxtries_per_branch = 40000,
  save_after_every_try = TRUE,
  savedir = getwd(),
  seedval = sample.int(1e9, 1),   
  wait_before_save = 0.01,
  master_nodenum_toPrint = 0
)

# ---- 5) Merge old + new ----
RES_clado_events_tables <- c(old_clado, BSM_output_new$RES_clado_events_tables)
RES_ana_events_tables   <- c(old_ana,   BSM_output_new$RES_ana_events_tables)

# ---- 6) Save merged 150-rep results ----
save(RES_clado_events_tables, file = paste0(model_name, "_RES_clado_events_tables_150.Rdata"))
save(RES_ana_events_tables,   file = paste0(model_name, "_RES_ana_events_tables_150.Rdata"))

BSM_output_150 <- list(
  RES_clado_events_tables = RES_clado_events_tables,
  RES_ana_events_tables   = RES_ana_events_tables
)
save(BSM_output_150, file = paste0(model_name, "_BSM_output_full_150.Rdata"))

# ---- 7) Quick verification ----
n_old <- length(old_clado)
n_new <- length(BSM_output_new$RES_clado_events_tables)
n_tot <- length(RES_clado_events_tables)

message(sprintf("Verification: old=%d, new=%d, total=%d (expected 150)", 
                n_old, n_new, n_tot))
stopifnot(n_tot == (n_old + n_new))
stopifnot(length(RES_ana_events_tables) == 
          (length(old_ana) + length(BSM_output_new$RES_ana_events_tables)))

message("\nSUCCESS: DEC+J BSM extended to 150 total replicates. Files written:\n",
        sprintf(" - %s_RES_clado_events_tables_150.Rdata\n", model_name),
        sprintf(" - %s_RES_ana_events_tables_150.Rdata\n", model_name),
        sprintf(" - %s_BSM_output_full_150.Rdata\n", model_name))

