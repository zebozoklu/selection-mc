# Project config
# Purpose: centralize global choices and defaults (seed, parameter defaults, paths)
# Inputs: none (or user edits)
# Outputs: a list/config object used by other scripts
# Sanity checks: paths exist; parameters have correct dimensions

make_config <- function() {
  list(
    # True outcome parameters (no intercept here; we can treat intercept separately)
    beta = c(1.0, -0.5),
    
    # Selection slopes for Z_noInt (currently same dimension as X)
    gamma_slopes = c(0.8, -0.4),
    
    # Which slope index is the "target" for sign reliability (1 means beta[1])
    target_beta_index = 1,
    
    # Monte Carlo reps
    R_debug = 50,
    R_full  = 1000,
    
    # Calibration sample size for Z
    n_calib = 200000,
    
    # Random seeds
    seed_calib = 1,
    seed_mc = 123
  )
}
