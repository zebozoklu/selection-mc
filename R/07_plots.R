# Plots
# Purpose: generate figures for bias/RMSE/sign-fail/coverage across scenarios
# Inputs: tidy summary table (output/results/mc_summary.rds)
# Outputs: saved figures in output/figs
# Sanity checks: consistent facet labels; readable axes; deterministic filenames

make_dirs_if_needed <- function() {
  if (!dir.exists("output")) dir.create("output")
  if (!dir.exists("output/figs")) dir.create("output/figs")
}

library(ggplot2)

make_plots <- function(summary_path = "output/results/mc_summary.rds",
                       drop_all_fail = TRUE,
                       warn_fail_rate_above = 0.20) {
  make_dirs_if_needed()
  
  if (!file.exists(summary_path)) {
    stop("Summary file not found: ", summary_path,
         "\nRun scripts/run_full_mc.R first.")
  }
  
  tab <- readRDS(summary_path)
  
  # focus on beta1 (x1)
  tab_beta1 <- subset(tab, coef_name == "x1")
  
  # drop scenarios where estimator never succeeded
  if (drop_all_fail && "n_used" %in% names(tab_beta1)) {
    tab_beta1 <- subset(tab_beta1, is.na(n_used) | n_used > 0)
  }
  
  # warn about high failure rates (informative for debugging/reporting)
  if ("fail_rate" %in% names(tab_beta1)) {
    worst <- aggregate(fail_rate ~ estimator, data = tab_beta1, FUN = max)
    worst <- worst[order(-worst$fail_rate), ]
    if (nrow(worst) > 0 && worst$fail_rate[1] >= warn_fail_rate_above) {
      message("Warning: high estimator failure rates detected (max by estimator):")
      print(worst)
    }
  }
  
  # nicer labels
  tab_beta1$rho_f   <- factor(tab_beta1$rho, levels = sort(unique(tab_beta1$rho)))
  tab_beta1$p_label <- paste0("p = ", tab_beta1$p_select)
  
  # show achieved selection rate if present (optional)
  if ("p_hat_mean" %in% names(tab_beta1)) {
    tab_beta1$p_label <- paste0(tab_beta1$p_label,
                                "\n(p̂ ≈ ", sprintf("%.2f", tab_beta1$p_hat_mean), ")")
  }
  
  tab_beta1$err_lab <- factor(tab_beta1$err_family,
                              levels = c("normal", "t3", "t5", "logistic"))
  
  tab_beta1$est_lab <- factor(
    tab_beta1$estimator,
    levels = c("selected_ols",
               "zero_impute",
               "heckman_probit",
               "heckman_logit",
               "heckman_lpm"),
    labels = c("Selected OLS",
               "Zero-impute",
               "Heckman (probit)",
               "Heckman (logit)",
               "Heckman (LPM)")
  )
  
  ## 1) Sign-failure plot for β1
  p_sign <- ggplot(tab_beta1,
                   aes(x = rho_f, y = sign_fail,
                       colour = est_lab, group = est_lab)) +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_point(position = position_dodge(width = 0.2)) +
    facet_grid(err_lab ~ p_label) +
    labs(
      x = expression(rho),
      y = "Sign failure rate for β[1]",
      colour = "Estimator",
      title = "Sign failure of β₁ across ρ, error family, and selection rate"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal()
  
  ggsave("output/figs/sign_fail_beta1_by_rho.png",
         p_sign, width = 9, height = 6, dpi = 300)
  
  ## 2) Bias plot for β1
  p_bias <- ggplot(tab_beta1,
                   aes(x = rho_f, y = bias,
                       colour = est_lab, group = est_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_point(position = position_dodge(width = 0.2)) +
    facet_grid(err_lab ~ p_label) +
    labs(
      x = expression(rho),
      y = "Bias of β[1]",
      colour = "Estimator",
      title = "Bias of β₁ across ρ, error family, and selection rate"
    ) +
    theme_minimal()
  
  ggsave("output/figs/bias_beta1_by_rho.png",
         p_bias, width = 9, height = 6, dpi = 300)
  
  ## 3) Coverage plot for β1
  p_cov <- ggplot(tab_beta1,
                  aes(x = rho_f, y = coverage,
                      colour = est_lab, group = est_lab)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_point(position = position_dodge(width = 0.2)) +
    facet_grid(err_lab ~ p_label) +
    labs(
      title  = "Coverage of 95% CIs for β₁ across ρ, error family, and selection rate",
      x      = expression(rho),
      y      = "Coverage probability",
      colour = "Estimator"
    ) +
    theme_minimal()
  
  ggsave("output/figs/coverage_beta1_by_rho.png",
         p_cov, width = 9, height = 6, dpi = 300)
  
  ## 4) (Optional) Failure-rate plot for β1 (same value across coefs, but plotted on x1 rows)
  p_fail <- NULL
  if ("fail_rate" %in% names(tab_beta1)) {
    p_fail <- ggplot(tab_beta1,
                     aes(x = rho_f, y = fail_rate,
                         colour = est_lab, group = est_lab)) +
      geom_line(position = position_dodge(width = 0.2)) +
      geom_point(position = position_dodge(width = 0.2)) +
      facet_grid(err_lab ~ p_label) +
      labs(
        title  = "Estimator failure rate across scenarios",
        x      = expression(rho),
        y      = "Failure rate",
        colour = "Estimator"
      ) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_minimal()
    
    ggsave("output/figs/fail_rate_by_rho.png",
           p_fail, width = 9, height = 6, dpi = 300)
  }
  
  invisible(list(
    sign_plot     = p_sign,
    bias_plot     = p_bias,
    coverage_plot = p_cov,
    fail_plot     = p_fail
  ))
}

if (sys.nframe() == 0L) {
  make_plots()
}



