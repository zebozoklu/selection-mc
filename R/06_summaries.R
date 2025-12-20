# R/06_summaries.R
# Performance summaries
# Purpose: compute bias, RMSE, sign-reversal frequency, coverage, and MCSEs
# Inputs: matrices of estimates + SEs + true parameter values
# Outputs: tidy summary table per scenario x estimator
# Notes: report Monte Carlo standard errors (MCSE)
#   - For probabilities (coverage, sign-fail): sqrt(p_hat*(1-p_hat)/R_eff)
#   - For Monte Carlo means: sd_emp / sqrt(R_eff)
#
# Key robustness improvements:
#   (i) MCSE denominators use effective sample sizes (non-missing reps)
#   (ii) coverage uses reps with finite beta_hat and finite se_hat
#   (iii) table includes beta_true (handy for density plots)

summarize_mc_estimator <- function(beta_hat_mat,
                                   se_hat_mat = NULL,
                                   beta_true,
                                   target_index = 1,
                                   zcrit = 1.96) {
  if (!is.matrix(beta_hat_mat)) stop("beta_hat_mat must be a matrix.")
  if (ncol(beta_hat_mat) != length(beta_true)) {
    stop("ncol(beta_hat_mat) must equal length(beta_true).")
  }
  
  R_total <- nrow(beta_hat_mat)
  K <- ncol(beta_hat_mat)
  
  # Coefficient names: ensure stable names exist (important for later plotting)
  coef_names <- colnames(beta_hat_mat)
  if (is.null(coef_names) || any(is.na(coef_names)) || any(coef_names == "")) {
    if (!is.null(names(beta_true)) && length(names(beta_true)) == K) {
      coef_names <- names(beta_true)
    } else {
      coef_names <- paste0("coef", seq_len(K))
    }
    colnames(beta_hat_mat) <- coef_names
    if (!is.null(se_hat_mat) && is.matrix(se_hat_mat)) colnames(se_hat_mat) <- coef_names
  }
  
  # Align true values with coefficient names
  if (is.null(names(beta_true))) names(beta_true) <- coef_names
  beta_true <- beta_true[coef_names]
  
  # Clean: treat non-finite as missing
  beta_use <- beta_hat_mat
  beta_use[!is.finite(beta_use)] <- NA_real_
  
  if (!is.null(se_hat_mat)) {
    if (!is.matrix(se_hat_mat)) stop("se_hat_mat must be a matrix when provided.")
    if (!all(dim(se_hat_mat) == dim(beta_hat_mat))) {
      stop("se_hat_mat must have same dimensions as beta_hat_mat.")
    }
    se_use <- se_hat_mat
    se_use[!is.finite(se_use)] <- NA_real_
    # (Optional) guard against negative SEs
    se_use[se_use < 0] <- NA_real_
  } else {
    se_use <- NULL
  }
  
  # Effective counts for beta (per coefficient)
  n_ok_beta <- colSums(!is.na(beta_use))
  
  # Monte Carlo mean and bias
  means <- colMeans(beta_use, na.rm = TRUE)
  bias  <- means - beta_true
  
  # RMSE
  diffs <- sweep(beta_use, 2, beta_true, FUN = "-")
  rmse  <- sqrt(colMeans(diffs^2, na.rm = TRUE))
  
  # Empirical SD of beta_hat
  sd_emp <- apply(beta_use, 2, sd, na.rm = TRUE)
  
  # MCSE of Monte Carlo mean: sd_emp / sqrt(R_eff)
  mcse_mean <- rep(NA_real_, K)
  ok_for_mcse <- n_ok_beta >= 2
  mcse_mean[ok_for_mcse] <- sd_emp[ok_for_mcse] / sqrt(n_ok_beta[ok_for_mcse])
  
  # ---- sign failure for the target coefficient (single scalar) + MCSE ----
  if (target_index < 1 || target_index > K) stop("target_index out of bounds.")
  
  target_coef_name <- coef_names[target_index]
  if (is.null(target_coef_name) || is.na(target_coef_name) || target_coef_name == "") {
    target_coef_name <- paste0("coef", target_index)
  }
  
  true_targ <- beta_true[target_index]
  targ <- beta_use[, target_index]
  n_ok_sign <- sum(!is.na(targ))
  
  if (!is.finite(true_targ) || true_targ == 0 || n_ok_sign == 0) {
    sign_fail <- NA_real_
    sign_fail_mcse <- NA_real_
  } else {
    # opposite sign (0 treated as not-a-sign-reversal)
    sign_fail <- mean(targ * true_targ < 0, na.rm = TRUE)
    sign_fail_mcse <- if (n_ok_sign > 0) sqrt(sign_fail * (1 - sign_fail) / n_ok_sign) else NA_real_
  }
  
  # ---- Coverage + SE diagnostics (if SE matrix provided) ----
  coverage <- rep(NA_real_, K)
  coverage_mcse <- rep(NA_real_, K)
  mean_se <- rep(NA_real_, K)
  n_ok_cov <- rep(0L, K)
  
  if (!is.null(se_use)) {
    mean_se <- colMeans(se_use, na.rm = TRUE)
    
    for (j in seq_len(K)) {
      valid <- !is.na(beta_use[, j]) & !is.na(se_use[, j])
      n_ok_cov[j] <- sum(valid)
      
      if (n_ok_cov[j] > 0) {
        lo <- beta_use[valid, j] - zcrit * se_use[valid, j]
        hi <- beta_use[valid, j] + zcrit * se_use[valid, j]
        
        cov_j <- (lo <= beta_true[j]) & (beta_true[j] <= hi)
        coverage[j] <- mean(cov_j)
        coverage_mcse[j] <- sqrt(coverage[j] * (1 - coverage[j]) / n_ok_cov[j])
      }
    }
  } else {
    n_ok_cov[] <- NA_integer_
  }
  
  # Ensure names
  names(means) <- coef_names
  names(bias) <- coef_names
  names(rmse) <- coef_names
  names(sd_emp) <- coef_names
  names(mcse_mean) <- coef_names
  names(coverage) <- coef_names
  names(coverage_mcse) <- coef_names
  names(mean_se) <- coef_names
  names(n_ok_beta) <- coef_names
  names(n_ok_cov) <- coef_names
  
  list(
    beta_true = beta_true,
    
    mean = means,
    bias = bias,
    rmse = rmse,
    
    sd_emp = sd_emp,
    mcse_mean = mcse_mean,
    
    # effective replication counts
    R_total = R_total,
    n_ok_beta = n_ok_beta,
    n_ok_cov  = n_ok_cov,
    
    sign_fail = sign_fail,                 # scalar
    sign_fail_mcse = sign_fail_mcse,       # scalar
    n_ok_sign = n_ok_sign,                 # scalar
    target_coef_name = target_coef_name,   # where sign_fail belongs
    
    coverage = coverage,
    coverage_mcse = coverage_mcse,
    mean_se = mean_se
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
    if (is.null(coef_names)) coef_names <- paste0("coef", seq_along(s$mean))
    
    # Put sign_fail only on the target coefficient row
    sign_fail_vec <- rep(NA_real_, length(coef_names))
    sign_fail_mcse_vec <- rep(NA_real_, length(coef_names))
    n_ok_sign_vec <- rep(NA_integer_, length(coef_names))
    
    if (!is.null(s$target_coef_name)) {
      hit <- which(coef_names == s$target_coef_name)
      if (length(hit) == 1L) {
        sign_fail_vec[hit] <- s$sign_fail
        sign_fail_mcse_vec[hit] <- s$sign_fail_mcse
        n_ok_sign_vec[hit] <- s$n_ok_sign
      }
    }
    
    tmp <- data.frame(
      scenario_id = scen$scenario_id,
      n           = scen$n,
      p_select    = scen$p_select,
      rho         = scen$rho,
      err_family  = as.character(scen$err_family),
      estimator   = est_name,
      coef_name   = coef_names,
      
      beta_true   = as.numeric(s$beta_true),
      
      mean        = as.numeric(s$mean),
      bias        = as.numeric(s$bias),
      rmse        = as.numeric(s$rmse),
      
      sd_emp      = as.numeric(s$sd_emp),
      mcse_mean   = as.numeric(s$mcse_mean),
      
      coverage       = as.numeric(s$coverage),
      coverage_mcse  = as.numeric(s$coverage_mcse),
      
      mean_se     = as.numeric(s$mean_se),
      
      # effective R info (helps interpret failures/nonconvergence)
      R_total     = as.integer(s$R_total),
      n_ok_beta   = as.integer(s$n_ok_beta),
      n_ok_cov    = as.integer(s$n_ok_cov),
      
      # sign_fail only for target coef; NA otherwise
      sign_fail      = sign_fail_vec,
      sign_fail_mcse = sign_fail_mcse_vec,
      n_ok_sign      = n_ok_sign_vec,
      
      stringsAsFactors = FALSE
    )
    
    rows[[est_name]] <- tmp
  }
  
  do.call(rbind, rows)
}


