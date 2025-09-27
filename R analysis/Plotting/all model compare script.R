# Clear and load pkgs
rm(list = ls())
library(BioGeoBEARS)
library(tidyverse)
library(gt)

# -----------------------------
# 1) Files & settings
# -----------------------------
files <- c(
  "DEC"            = "res_DEC.Rdata",
  "DEC+J"          = "res_DECJ.Rdata",
  "BAYAREALIKE"    = "BAYAREALIKE_results.Rdata",
  "BAYAREALIKE+J"  = "BAYAREALIKEJ_results.Rdata",
  "DIVALIKE"       = "DIVALIKE_results.Rdata",
  "DIVALIKE+J"     = "DIVALIKEJ_results.Rdata"
)

# Number of tips used in your BioGeoBEARS run
n_tips <- 351  # <-- change if needed

# -----------------------------
# 2) Helpers
# -----------------------------
fmt_param <- function(x) {
  if (is.na(x)) return("NA")
  if (x == 0) return("0")
  if (abs(x) < 1e-4) {
    formatC(x, format = "e", digits = 2)   # e.g., 2.61e-06
  } else if (abs(x) < 1) {
    as.character(signif(x, 3))             # 3 sig figs for small values
  } else {
    as.character(round(x, 3))              # 3 dp for larger values
  }
}

# Extract estimate and annotate if at max bound
pull_param <- function(pt, pname) {
  if (!pname %in% rownames(pt)) return(NA_character_)
  est <- suppressWarnings(as.numeric(pt[pname, "est"]))
  mx  <- suppressWarnings(as.numeric(pt[pname, "max"]))
  if (is.na(est)) return(NA_character_)
  out <- fmt_param(est)
  # mark “(max)” if near upper bound
  if (!is.na(mx) && !is.na(est) && (est >= (mx - 1e-12))) {
    out <- paste0(out, " (max)")
  }
  out
}

safe_num <- function(x) { suppressWarnings(as.numeric(x)) }

# -----------------------------
# 3) Load, extract, compute
# -----------------------------
rows <- list()

for (m in names(files)) {
  f <- files[[m]]
  if (!file.exists(f)) {
    warning(paste("File not found:", f))
    rows[[m]] <- tibble(
      Model = m, LnL = NA_real_, k = NA_integer_, AIC = NA_real_, AICc = NA_real_,
      d = NA_character_, e = NA_character_, j = NA_character_
    )
    next
  }
  env <- new.env()
  load(f, envir = env)
  
  # Guess the object name
  obj_name <- intersect(ls(env), c("res", "resDEC", "resBAYAREALIKE", "resBAYAREALIKEJ",
                                   "resDIVALIKE", "resDIVALIKEJ"))
  if (length(obj_name) == 0) {
    # fallback: take the first object in the file
    obj_name <- ls(env)[1]
  }
  res_obj <- env[[obj_name]]
  
  # LnL
  lnL <- NA_real_
  if (!is.null(res_obj$total_loglikelihood)) {
    lnL <- safe_num(res_obj$total_loglikelihood)
  }
  
  # params table
  pt <- tryCatch(res_obj$outputs@params_table, error = function(e) NULL)
  
  # number of free parameters k (robust across models)
  k <- NA_integer_
  if (!is.null(pt) && "type" %in% colnames(pt)) {
    k <- sum(pt[, "type", drop = TRUE] == "free", na.rm = TRUE)
  } else {
    # conservative fallback: DEC=2, +J=3
    k <- ifelse(grepl("\\+J$", m), 3L, 2L)
  }
  
  # AIC & AICc
  AIC  <- if (!is.na(lnL) && !is.na(k)) 2*k - 2*lnL else NA_real_
  AICc <- if (!is.na(AIC) && !is.na(k) && !is.na(n_tips) && (n_tips - k - 1) > 0) {
    AIC + (2 * k * (k + 1)) / (n_tips - k - 1)
  } else NA_real_
  
  # Parameters d/e/j (as strings, nicely formatted)
  d_str <- if (!is.null(pt)) pull_param(pt, "d") else NA_character_
  e_str <- if (!is.null(pt)) pull_param(pt, "e") else NA_character_
  j_str <- if (!is.null(pt)) pull_param(pt, "j") else NA_character_
  
  rows[[m]] <- tibble(
    Model = m, LnL = lnL, k = k, AIC = AIC, AICc = AICc,
    d = d_str, e = e_str, j = j_str
  )
}

model_df <- bind_rows(rows) |>
  arrange(AICc)

# winner for highlighting
winner <- model_df$Model[which.min(model_df$AICc)]

# -----------------------------
# 4) Pretty table (gt)
# -----------------------------
model_df %>%
  select(Model, LnL, AICc, d, e, j) %>%
  gt() %>%
  tab_header(
    title = "Model Comparison: BioGeoBEARS",
    subtitle = "Sorted by AICc (lower is better)"
  ) %>%
  fmt_number(columns = c(LnL, AICc), decimals = 2) %>%
  cols_label(
    Model = "Model",
    LnL   = "Log-Likelihood",
    AICc  = "AICc",
    d     = "d (dispersal)",
    e     = "e (extinction)",
    j     = "j (founder-event)"
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#D9EAD3"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(rows = Model == winner)
  ) %>%
  tab_footnote(
    footnote = paste0("n (tips) = ", n_tips, 
                      "; k counted as number of free parameters in BioGeoBEARS params table.")
  ) %>%
  gtsave("Model_Comparison_AllModels.png")
