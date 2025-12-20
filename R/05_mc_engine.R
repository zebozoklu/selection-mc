# Monte Carlo engine
# Purpose: repeat simulation + estimation R times, for each scenario
# Inputs: scenario table, number of reps R, estimator functions
# Outputs:
#   - raw_beta: estimates per rep (R x K) per estimator
#   - raw_se:   robust SEs per rep (R x K) per estimator
#   - aggregated summaries per scenario/estimator
# Sanity checks: reproducible with seed; progress logging; saves intermediate results

run_mc_one_scenario <- function(scen, cfg, estimators,
                                R = cfg$R_debug,
                                seed = cfg$seed_mc) {
  
  set.seed(seed)
  
  p <- length(cfg$beta)       # number of X slopes
  K <- p + 1                  # intercept + slopes
  
  # prepare containers
  raw_beta <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K)
  })
  raw_se <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K)
  })
  
  # optional diagnostic
  sel_rate <- rep(NA_real_, R)
  
  # main loop
  for (r in seq_len(R)) {
    dat <- simulate_one_dataset(scen, cfg)
    
    sel_rate[r] <- mean(dat$S)
    
    for (nm in names(estimators)) {
      out <- estimators[[nm]](dat)
      
      if (!is.list(out) || is.null(out$beta) || is.null(out$se)) {
        stop("Estimator ", nm, " must return list(beta=..., se=...).")
      }
      
      b_hat <- out$beta
      se_hat <- out$se
      
      if (length(b_hat) != K || length(se_hat) != K) {
        stop("Estimator ", nm, " returned beta/se length ",
             length(b_hat), "/", length(se_hat),
             " but expected ", K, "/", K, ".")
      }
      
      # On first iteration, set column names from b_hat
      if (r == 1L) {
        colnames(raw_beta[[nm]]) <- names(b_hat)
        colnames(raw_se[[nm]]) <- names(se_hat)
      }
      
      raw_beta[[nm]][r, ] <- b_hat
      raw_se[[nm]][r, ] <- se_hat
    }
  }
  
  # true coefficients (we assume true intercept = 0)
  beta_true <- c(0, cfg$beta)
  
  # summarize (we'll update summarize_mc_estimator to optionally accept se draws)
  summary <- lapply(names(raw_beta), function(nm) {
    summarize_mc_estimator(
      beta_hat_mat = raw_beta[[nm]],
      se_hat_mat   = raw_se[[nm]],
      beta_true    = beta_true,
      target_index = 1 + cfg$target_beta_index
    )
  })
  names(summary) <- names(raw_beta)
  
  list(
    raw_beta = raw_beta,
    raw_se   = raw_se,
    sel_rate = sel_rate,
    summary  = summary,
    scen     = scen
  )
}



