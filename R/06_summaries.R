# Performance summaries
# Purpose: compute bias, RMSE, sign-reversal frequency (and maybe coverage later)
# Inputs: matrix/array of estimates + true parameter values
# Outputs: tidy summary table per scenario x estimator
# Sanity checks: metrics computed on correct coefficient index; handles NA gracefully

summarize_mc_estimator <- function(est_mat, beta_true, target_index = 1) {
  # est_mat: R x K matrix of estimates (rows = replications, cols = coefficients)
  # beta_true: length-K vector of true coefficients (intercept + slopes)
  # target_index: index of the slope in beta_true we care about for sign failures
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
  
  list(
    mean = means,
    bias = bias,
    rmse = rmse,
    sign_fail = sign_fail
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
      stringsAsFactors = FALSE
    )
    
    rows[[est_name]] <- tmp
  }
  
  do.call(rbind, rows)
}
