#######################################################
# BioGeoBEARS: Run BSM on DEC+J model and save outputs
#######################################################

library(BioGeoBEARS)

# 1. Load DEC+J results object from HPC run
load("DECJ_results.Rdata")  # This must contain `resDECJ`
res <- resDECj
model_name <- "DECJ"

# 2. Get inputs for stochastic mapping
stochastic_mapping_inputs_list <- get_inputs_for_stochastic_mapping(res = res)
save(stochastic_mapping_inputs_list, file = paste0(model_name, "_BSM_inputs.Rdata"))

# 3. Run the BSM (50 maps)
BSM_output <- runBSM(
  res,
  stochastic_mapping_inputs_list = stochastic_mapping_inputs_list,
  maxnum_maps_to_try = 100,
  nummaps_goal = 50,
  maxtries_per_branch = 40000,
  save_after_every_try = TRUE,
  savedir = getwd(),
  seedval = 12345,
  wait_before_save = 0.01,
  master_nodenum_toPrint = 0
)

# 4. Extract and save outputs for later plotting in RStudio GUI
RES_clado_events_tables <- BSM_output$RES_clado_events_tables
RES_ana_events_tables <- BSM_output$RES_ana_events_tables

save(RES_clado_events_tables, file = paste0(model_name, "_RES_clado_events_tables.Rdata"))
save(RES_ana_events_tables, file = paste0(model_name, "_RES_ana_events_tables.Rdata"))

# 5. Save the BSM_output object itself (optional, contains everything)
save(BSM_output, file = paste0(model_name, "_BSM_output_full.Rdata"))

cat("\nBSM complete for", model_name, "– outputs saved for local plotting.\n")

