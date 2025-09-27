#######################################################
# Run BAYAREALIKE+J BSMs to a total of 150 replicates
#######################################################

library(BioGeoBEARS)

model_name <- "BAYAREALIKEJ"
results_rdata <- "BAYAREALIKEJ_results.Rdata"  # must contain: resBAYAREALIKEj

# ---- 1) Load fitted results ----
load(results_rdata)        # loads resBAYAREALIKEj
res <- resBAYAREALIKEj

# ---- 2) Load any existing BSMs (if you’ve run some already) ----
old_clado <- old_ana <- NULL

if (file.exists(paste0(model_name, "_BSM_output_full.Rdata"))) {
  load(paste0(model_name, "_BSM_output_full.Rdata"))
  if (exists("BSM_output")) {
    old_clado <- BSM_output$RES_clado_events_tables
    old_ana   <- BSM_output$RES_ana_events_tables
  }
}
if (is.null(old_clado) && file.exists(paste0(model_name, "_RES_clado_events_tables.Rdata"))) {
  load(paste0(model_name, "_RES_clado_events_tables.Rdata"))
  old_clado <- RES_clado_events_tables
}
if (is.null(old_ana) && file.exists(paste0(model_name, "_RES_ana_events_tables.Rdata"))) {
  load(paste0(model_name, "_RES_ana_events_tables.Rdata"))
  old_ana <- RES_ana_events_tables
}

n_old <- if (is.null(old_clado)) 0 else length(old_clado)
message(sprintf("Found %d existing BSMs for %s.", n_old, model_name))

# ---- 3) Target and needed ----
target_total <- 150L
n_to_add <- max(0L, target_total - n_old)
if (n_to_add == 0L) {
  message("Already have >= target BSMs. Nothing to do.")
  quit(save="no")
}

# ---- 4) Prepare & run ----
sminputs <- get_inputs_for_stochastic_mapping(res = res)

set.seed(654321)
BSM_new <- runBSM(
  res,
  stochastic_mapping_inputs_list = sminputs,
  maxnum_maps_to_try = n_to_add,
  nummaps_goal = n_to_add,
  maxtries_per_branch = 40000,
  save_after_every_try = TRUE,
  savedir = getwd(),
  seedval = sample.int(1e9, 1),
  wait_before_save = 0.01,
  master_nodenum_toPrint = 0
)

# ---- 5) Merge + save ----
RES_clado_events_tables <- c(if (!is.null(old_clado)) old_clado else list(),
                             BSM_new$RES_clado_events_tables)
RES_ana_events_tables   <- c(if (!is.null(old_ana))   old_ana   else list(),
                             BSM_new$RES_ana_events_tables)

save(RES_clado_events_tables, file = paste0(model_name, "_RES_clado_events_tables_150.Rdata"))
save(RES_ana_events_tables,   file = paste0(model_name, "_RES_ana_events_tables_150.Rdata"))

BSM_output_150 <- list(
  RES_clado_events_tables = RES_clado_events_tables,
  RES_ana_events_tables   = RES_ana_events_tables
)
save(BSM_output_150, file = paste0(model_name, "_BSM_output_full_150.Rdata"))

n_tot <- length(RES_clado_events_tables)
message(sprintf("SUCCESS: %s BSMs now at %d total (target %d).",
                model_name, n_tot, target_total))

