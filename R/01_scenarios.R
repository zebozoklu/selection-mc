# Scenario grid
# Purpose: define the Monte Carlo design grid (n, selection rate p, rho, error family, sel_model)
# Inputs: user-chosen vectors of knob values
# Outputs: data.frame/tibble with one row per scenario + scenario_id
# Sanity checks: grid size is expected; no duplicated scenario_id

make_scenario_grid <- function(
    n_vec        = 500,                          # single n
    p_select_vec = c(0.30, 0.60),                # two selection rates
    rho_vec      = c(-0.6, 0.0, 0.6),            # negative, zero, positive selection on unobservables
    err_families = c("normal", "t3", "t5", "logistic"),
    sel_models   = c("probit", "logit", "lpm")   # true selection links
) {
  grid <- expand.grid(
    n         = n_vec,
    p_select  = p_select_vec,
    rho       = rho_vec,
    err_family = err_families,
    sel_model  = sel_models,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # attach error parameters (df for t-distributions; logistic has no extra parameter)
  grid$df <- NA_integer_
  grid$df[grid$err_family == "t3"] <- 3L
  grid$df[grid$err_family == "t5"] <- 5L
  
  # nice unique ID
  grid$scenario_id <- sprintf("S%03d", seq_len(nrow(grid)))
  
  # reorder columns for readability
  grid <- grid[, c("scenario_id",
                   "n",
                   "p_select",
                   "rho",
                   "err_family",
                   "sel_model",
                   "df")]
  
  grid
}


