# Monte Carlo engine
# Purpose: repeat simulation + estimation R times, for each scenario
# Inputs: scenario row, config, estimators, reps R, seed
# Outputs:
#   - raw draws (estimates per rep) per estimator
#   - raw SEs per rep per estimator
#   - failure stats per estimator
#   - scenario diagnostics (p_hat, uv_corr)
#   - aggregated summaries per scenario/estimator (computed on successful reps only)

run_mc_one_scenario <- function(scen, cfg, estimators,
                                R = cfg$R_debug,
                                seed = cfg$seed_mc) {
  
  set.seed(seed)
  
  p <- length(cfg$beta)  # number of X slopes
  K <- p + 1             # intercept + slopes
  
  # containers for betas and ses
  raw    <- lapply(estimators, function(f) matrix(NA_real_, nrow = R, ncol = K))
  raw_se <- lapply(estimators, function(f) matrix(NA_real_, nrow = R, ncol = K))
  
  # failure tracking: one logical vector per estimator (TRUE if failed at rep r)
  failed <- lapply(estimators, function(f) rep(FALSE, R))
  
  # scenario diagnostics (per rep)
  p_hat_rep   <- rep(NA_real_, R)
  uv_corr_rep <- rep(NA_real_, R)
  
  # main loop
  for (r in seq_len(R)) {
    dat <- simulate_one_dataset(scen, cfg)
    
    # collect DGP diagnostics (if present)
    if (!is.null(dat$p_hat))   p_hat_rep[r]   <- dat$p_hat
    if (!is.null(dat$uv_corr)) uv_corr_rep[r] <- dat$uv_corr
    
    for (nm in names(estimators)) {
      res <- tryCatch(
        estimators[[nm]](dat, cfg),
        error = function(e) e
      )
      
      if (inherits(res, "error")) {
        failed[[nm]][r] <- TRUE
        next
      }
      
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
      
      # On first successful fill, set column names if they aren't set
      if (is.null(colnames(raw[[nm]]))) {
        colnames(raw[[nm]])    <- names(beta_vec)
        colnames(raw_se[[nm]]) <- names(beta_vec)
      }
      
      raw[[nm]][r, ]    <- beta_vec
      raw_se[[nm]][r, ] <- se_vec
    }
  }
  
  # true coefficients (we assume true intercept = 0)
  beta_true <- c(0, cfg$beta)
  
  # scenario-level diagnostics
  diag <- list(
    p_hat_mean   = mean(p_hat_rep,   na.rm = TRUE),
    p_hat_sd     = sd(p_hat_rep,     na.rm = TRUE),
    uv_corr_mean = mean(uv_corr_rep, na.rm = TRUE),
    uv_corr_sd   = sd(uv_corr_rep,   na.rm = TRUE)
  )
  
  # summarize: one summary list per estimator (successful reps only)
  summary <- list()
  fail_stats <- data.frame(
    estimator = names(estimators),
    n_fail    = NA_integer_,
    fail_rate = NA_real_,
    n_used    = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  for (j in seq_along(estimators)) {
    nm <- names(estimators)[j]
    
    ok <- !failed[[nm]]
    n_fail <- sum(!ok)
    n_used <- sum(ok)
    
    fail_stats$n_fail[j]    <- n_fail
    fail_stats$fail_rate[j] <- n_fail / R
    fail_stats$n_used[j]    <- n_used
    
    if (n_used == 0L) {
      # everything failed: return NA summary but keep fail stats
      summary[[nm]] <- list(
        mean      = rep(NA_real_, K),
        bias      = rep(NA_real_, K),
        rmse      = rep(NA_real_, K),
        sign_fail = NA_real_,
        coverage  = NA_real_
      )
      next
    }
    
    summary[[nm]] <- summarize_mc_estimator(
      est_mat      = raw[[nm]][ok, , drop = FALSE],
      beta_true    = beta_true,
      target_index = 1 + cfg$target_beta_index,
      se_mat       = raw_se[[nm]][ok, , drop = FALSE]
    )
    
    # preserve coefficient names if available
    if (!is.null(colnames(raw[[nm]]))) {
      names(summary[[nm]]$mean) <- colnames(raw[[nm]])
      names(summary[[nm]]$bias) <- colnames(raw[[nm]])
      names(summary[[nm]]$rmse) <- colnames(raw[[nm]])
    }
  }
  
  list(
    raw        = raw,
    raw_se     = raw_se,
    failed     = failed,
    fail_stats = fail_stats,
    diag       = diag,
    summary    = summary,
    scen       = scen
  )
}



