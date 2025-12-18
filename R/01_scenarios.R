# Scenario grid
# Purpose: define the Monte Carlo design grid (n, selection rate p, rho, error family)
# Inputs: user-chosen vectors of knob values
# Outputs: data.frame/tibble with one row per scenario + scenario_id
# Sanity checks: grid size is expected; no duplicated scenario_id


make_scenario_grid <- function(
    n_vec = c(500, 2000),
    p_select_vec = c(0.30, 0.60),
    rho_vec = c(0.0, 0.3, 0.6),
    err_families = c("normal", "t3", "t5", "mixture")
) {
  grid <- expand.grid(
    n = n_vec,
    p_select = p_select_vec,
    rho = rho_vec,
    err_family = err_families,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # attach error parameters (keep as explicit columns for simplicity)
  grid$df <- NA_integer_
  grid$mix_pi <- NA_real_
  grid$mix_shift_u <- NA_real_
  
  grid$df[grid$err_family == "t3"] <- 3L
  grid$df[grid$err_family == "t5"] <- 5L
  
  grid$mix_pi[grid$err_family == "mixture"] <- 0.10
  grid$mix_shift_u[grid$err_family == "mixture"] <- 2.0
  
  # nice unique ID
  grid$scenario_id <- sprintf("S%03d", seq_len(nrow(grid)))
  
  # reorder columns
  grid <- grid[, c("scenario_id","n","p_select","rho","err_family","df","mix_pi","mix_shift_u")]
  
  grid
}
