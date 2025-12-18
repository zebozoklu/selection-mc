# Project config
# Purpose: centralize global choices and defaults (seed, parameter defaults, paths)
# Inputs: none (or user edits)
# Outputs: a list/config object used by other scripts
# Sanity checks: paths exist; parameters have correct dimensions

`%||%` <- function(a, b) if (!is.null(a)) a else b

make_config <- function() {
  list(
    ## -----------------------
    ## True parameters
    ## -----------------------
    # Outcome slopes (no intercept; true intercept assumed 0)
    beta = c(0.05, -0.025),
    
    # Selection slopes on X (Z_noInt). If later you add an exclusion C,
    # you'll expand this structure.
    gamma_slopes = c(0.8, -0.4),
    
    # Target coefficient index among beta slopes (1 means beta[1] = x1)
    target_beta_index = 1,
    
    ## -----------------------
    ## Covariate distribution
    ## -----------------------
    # X ~ N(0, SigmaX) where SigmaX has AR correlation
    SigmaX_ar = 0.3,
    
    ## -----------------------
    ## Inference / SE settings
    ## -----------------------
    # "HC1" recommended (robust sandwich). Optionally "homosked".
    se_type = "HC1",
    
    # numerical stability for probabilities
    eps_prob = 1e-8,
    
    ## -----------------------
    ## Monte Carlo settings
    ## -----------------------
    R_debug = 50,
    R_full  = 1000,
    
    # Calibration sample size
    n_calib = 200000,
    
    ## -----------------------
    ## Seeds
    ## -----------------------
    seed_calib = 1,
    seed_mc    = 123
  )
}




