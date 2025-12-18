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
  # keep coefficient names if present
  if (!is.null(colnames(est_mat))) {
    names(means) <- colnames(est_mat)
  }
  
  bias <- means - beta_true
  
  # center each column by true value then square
  diffs <- sweep(est_mat, 2, beta_true, FUN = "-")
  rmse  <- sqrt(colMeans(diffs^2))
  
  names(bias) <- names(means)
  names(rmse) <- names(means)
  
  # sign failure for the target coefficient
  sign_fail <- mean(sign(est_mat[, target_index]) != sign(beta_true[target_index]))
  
  # coverage for target coefficient, if SEs available
  coverage <- NA_real_
  if (!is.null(se_mat)) {
    if (!is.matrix(se_mat) ||
        any(dim(se_mat) != dim(est_mat))) {
      stop("se_mat must be a matrix with same dim as est_mat.")
    }
    zcrit <- qnorm(1 - alpha / 2)
    est_t <- est_mat[, target_index]
    se_t  <- se_mat[, target_index]
    
    lower <- est_t - zcrit * se_t
    upper <- est_t + zcrit * se_t
    
    inside <- (beta_true[target_index] >= lower) &
      (beta_true[target_index] <= upper)
    coverage <- mean(inside)
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
  scen <- mc_res$scen
  summaries <- mc_res$summary
  
  rows <- list()
  
  for (est_name in names(summaries)) {
    s <- summaries[[est_name]]
    
    coef_names <- names(s$mean)
    if (is.null(coef_names)) {
      coef_names <- paste0("coef", seq_along(s$mean))
    }
    
    tmp <- data.frame(
      scenario_id = scen$scenario_id,
      n           = scen$n,
      p_select    = scen$p_select,
      rho         = scen$rho,
      err_family  = as.character(scen$err_family),
      estimator   = est_name,
      coef_name   = coef_names,
      mean        = as.numeric(s$mean),
      bias        = as.numeric(s$bias),
      rmse        = as.numeric(s$rmse),
      sign_fail   = s$sign_fail,
      coverage    = s$coverage,   # same value for all coef rows
      stringsAsFactors = FALSE
    )
    
    rows[[est_name]] <- tmp
  }
  
  do.call(rbind, rows)
}



