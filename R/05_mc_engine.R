# Monte Carlo engine
# Purpose: repeat simulation + estimation R times, for each scenario
# Inputs: scenario table, number of reps R, estimator functions
# Outputs:
#   - raw draws (estimates per rep) OR
#   - aggregated summaries per scenario/estimator
# Sanity checks: reproducible with seed; progress logging; saves intermediate results

run_mc_one_scenario <- function(scen, cfg, estimators,
                                R = cfg$R_debug,
                                seed = cfg$seed_mc) {
  
  set.seed(seed)
  
  p <- length(cfg$beta)       # number of X slopes
  K <- p + 1                  # intercept + slopes
  
  # prepare containers
  raw <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K)
  })
  
  # main loop
  for (r in seq_len(R)) {
    dat <- simulate_one_dataset(scen, cfg)
    
    for (nm in names(estimators)) {
      b_hat <- estimators[[nm]](dat)
      
      if (length(b_hat) != K) {
        stop("Estimator ", nm, " returned length ", length(b_hat),
             " but expected ", K, ".")
      }
      
      raw[[nm]][r, ] <- b_hat
    }
  }
  
  # true coefficients (we assume true intercept = 0)
  beta_true <- c(0, cfg$beta)
  
  # summarize
  summary <- lapply(raw, summarize_mc_estimator,
                    beta_true = beta_true,
                    target_index = 1 + cfg$target_beta_index)
  
  list(
    raw = raw,
    summary = summary,
    scen = scen
  )
}
