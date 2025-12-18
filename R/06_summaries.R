# Performance summaries
# Purpose: compute bias, RMSE, sign-reversal freq, and coverage
# Inputs: est_mat (R x K), beta_true (length K), optional se_mat (R x K)
# Outputs: list(mean, bias, rmse, sign_fail, coverage)
# coverage is for the target coefficient only (scalar)

summarize_mc_estimator <- function(est_mat, beta_true, target_index = 1,
                                   se_mat = NULL, alpha = 0.05) {
  if (!is.matrix(est_mat)) stop("est_mat must be a matrix.")
  if (ncol(est_mat) != length(beta_true)) {
    stop("ncol(est_mat) must equal length(beta_true).")
  }
  
  means <- colMeans(est_mat)
  if (!is.null(colnames(est_mat))) names(means) <- colnames(est_mat)
  
  bias <- means - beta_true
  
  diffs <- sweep(est_mat, 2, beta_true, FUN = "-")
  rmse  <- sqrt(colMeans(diffs^2))
  
  names(bias) <- names(means)
  names(rmse) <- names(means)
  
  # Sign failure (clean definition):
  # count strict sign reversals relative to the true sign
  bt <- beta_true[target_index]
  if (bt == 0) {
    sign_fail <- NA_real_
  } else {
    sign_fail <- mean(est_mat[, target_index] * bt < 0)
  }
  
  # Coverage for target coefficient, if SEs available
  coverage <- NA_real_
  if (!is.null(se_mat)) {
    if (!is.matrix(se_mat) || any(dim(se_mat) != dim(est_mat))) {
      stop("se_mat must be a matrix with same dim as est_mat.")
    }
    
    zcrit <- qnorm(1 - alpha / 2)
    est_t <- est_mat[, target_index]
    se_t  <- se_mat[, target_index]
    
    lower <- est_t - zcrit * se_t
    upper <- est_t + zcrit * se_t
    
    inside <- (beta_true[target_index] >= lower) &
      (beta_true[target_index] <= upper)
    coverage <- mean(inside, na.rm = TRUE)
  }
  
  list(
    mean      = means,
    bias      = bias,
    rmse      = rmse,
    sign_fail = sign_fail,
    coverage  = coverage
  )
}

# Turn one scenario's MC result into a tidy summary table
# mc_res: output of run_mc_one_scenario()
mc_result_to_table <- function(mc_res) {
  scen      <- mc_res$scen
  summaries <- mc_res$summary
  
  # scenario-level diagnostics (from MC engine)
  diag <- mc_res$diag %||% list()
  p_hat_mean   <- diag$p_hat_mean
  p_hat_sd     <- diag$p_hat_sd
  uv_corr_mean <- diag$uv_corr_mean
  uv_corr_sd   <- diag$uv_corr_sd
  
  # estimator-level failure stats
  fail_stats <- mc_res$fail_stats
  if (is.null(fail_stats) || !is.data.frame(fail_stats)) {
    fail_stats <- data.frame(
      estimator = names(summaries),
      n_fail    = NA_integer_,
      fail_rate = NA_real_,
      n_used    = NA_integer_,
      stringsAsFactors = FALSE
    )
  }
  
  rows <- list()
  
  for (est_name in names(summaries)) {
    s <- summaries[[est_name]]
    
    coef_names <- names(s$mean)
    if (is.null(coef_names)) coef_names <- paste0("coef", seq_along(s$mean))
    
    fs <- fail_stats[fail_stats$estimator == est_name, , drop = FALSE]
    if (nrow(fs) == 0L) {
      n_fail    <- NA_integer_
      fail_rate <- NA_real_
      n_used    <- NA_integer_
    } else {
      n_fail    <- as.integer(fs$n_fail[1])
      fail_rate <- as.numeric(fs$fail_rate[1])
      n_used    <- as.integer(fs$n_used[1])
    }
    
    tmp <- data.frame(
      scenario_id = scen$scenario_id,
      n           = scen$n,
      p_select    = scen$p_select,
      rho         = scen$rho,
      err_family  = as.character(scen$err_family),
      
      # NEW: diagnostics / hats
      p_hat_mean   = p_hat_mean,
      p_hat_sd     = p_hat_sd,
      uv_corr_mean = uv_corr_mean,
      uv_corr_sd   = uv_corr_sd,
      
      estimator   = est_name,
      n_used      = n_used,
      n_fail      = n_fail,
      fail_rate   = fail_rate,
      
      coef_name   = coef_names,
      mean        = as.numeric(s$mean),
      bias        = as.numeric(s$bias),
      rmse        = as.numeric(s$rmse),
      sign_fail   = s$sign_fail,
      coverage    = s$coverage,   # repeated across coef rows (target coef only)
      
      stringsAsFactors = FALSE
    )
    
    rows[[est_name]] <- tmp
  }
  
  do.call(rbind, rows)
}



