#### KISS: extend DEC+J BSMs to 500 total ---------------------------------

library(BioGeoBEARS)

model_name <- "DECJ"
target_total <- 500L

# 1) Load fitted results (must contain resDECj)
load("DECJ_results.Rdata"); res <- resDECj

# 2) Load existing 150 BSMs (prefer the bundled file you already have)
old_clado <- old_ana <- NULL
if (file.exists("DECJ_BSM_output_full_150.Rdata")) {
  load("DECJ_BSM_output_full_150.Rdata")  # creates BSM_output_150
  x <- if (exists("BSM_output_150")) BSM_output_150 else BSM_output
  old_clado <- x$RES_clado_events_tables
  old_ana   <- x$RES_ana_events_tables
} else {
  load("DECJ_RES_clado_events_tables_150.Rdata"); old_clado <- RES_clado_events_tables
  load("DECJ_RES_ana_events_tables_150.Rdata");   old_ana   <- RES_ana_events_tables
}

n_old <- length(old_clado)
message(sprintf("Existing DEC+J BSMs: %d", n_old))

# 3) How many more?
n_to_add <- max(0L, target_total - n_old)
if (n_to_add == 0L) { message("Nothing to add."); quit(save="no") }

# 4) Prepare inputs + run the extra maps
sminputs <- get_inputs_for_stochastic_mapping(res = res)

set.seed(54321)
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

# 5) Merge + save as *_500
RES_clado_events_tables <- c(old_clado, BSM_new$RES_clado_events_tables)
RES_ana_events_tables   <- c(old_ana,   BSM_new$RES_ana_events_tables)

stopifnot(length(RES_clado_events_tables) == target_total)
stopifnot(length(RES_ana_events_tables)   == target_total)

save(RES_clado_events_tables, file = paste0(model_name, "_RES_clado_events_tables_500.Rdata"))
save(RES_ana_events_tables,   file = paste0(model_name, "_RES_ana_events_tables_500.Rdata"))
BSM_output_500 <- list(RES_clado_events_tables=RES_clado_events_tables,
                       RES_ana_events_tables=RES_ana_events_tables)
save(BSM_output_500, file = paste0(model_name, "_BSM_output_full_500.Rdata"))

message("SUCCESS: DEC+J BSMs now at 500 total.")

