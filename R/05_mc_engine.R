# Monte Carlo engine
# Purpose: repeat simulation + estimation R times, for each scenario
# Inputs: scenario table, number of reps R, estimator functions
# Outputs:
#   - raw draws (estimates per rep)
#   - raw SEs per rep
#   - aggregated summaries per scenario/estimator
# Sanity checks: reproducible with seed; progress logging

run_mc_one_scenario <- function(scen, cfg, estimators,
                                R = cfg$R_debug,
                                seed = cfg$seed_mc) {
  
  set.seed(seed)
  
  p <- length(cfg$beta)       # number of X slopes
  K <- p + 1                  # intercept + slopes
  
  # containers for betas and ses
  raw    <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K)
  })
  raw_se <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K)
  })
  
  # main loop
  for (r in seq_len(R)) {
    dat <- simulate_one_dataset(scen, cfg)
    
    for (nm in names(estimators)) {
      res <- estimators[[nm]](dat)
      
      # allow either list(beta,se) or plain numeric (fallback)
      if (is.list(res) && !is.null(res$beta)) {
        beta_vec <- res$beta
        se_vec   <- res$se
      } else {
        beta_vec <- as.numeric(res)
        se_vec   <- rep(NA_real_, length(beta_vec))
      }
      
      if (length(beta_vec) != K) {
        stop("Estimator ", nm, " returned length ", length(beta_vec),
             " but expected ", K, ".")
      }
      
      # On first iteration, set column names from beta_vec
      if (r == 1L) {
        colnames(raw[[nm]])    <- names(beta_vec)
        colnames(raw_se[[nm]]) <- names(beta_vec)
      }
      
      raw[[nm]][r, ]    <- beta_vec
      raw_se[[nm]][r, ] <- se_vec
    }
  }
  
  # true coefficients (we assume true intercept = 0)
  beta_true <- c(0, cfg$beta)
  
  # summarize: one summary list per estimator
  summary <- list()
  for (nm in names(raw)) {
    summary[[nm]] <- summarize_mc_estimator(
      est_mat     = raw[[nm]],
      beta_true   = beta_true,
      target_index = 1 + cfg$target_beta_index,
      se_mat      = raw_se[[nm]]
    )
  }
  
  list(
    raw    = raw,
    raw_se = raw_se,
    summary = summary,
    scen    = scen
  )
}




