# Reporting script
# Goal:
#   - load saved results
#   - produce summary tables (bias/RMSE/sign-fail/coverage)
#   - produce and save plots to output/figs/

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

make_dirs_if_needed <- function() {
  if (!dir.exists("output")) dir.create("output")
  if (!dir.exists("output/figs")) dir.create("output/figs")
  if (!dir.exists("output/results")) dir.create("output/results", recursive = TRUE)
  if (!dir.exists("output/results/tables")) dir.create("output/results/tables", recursive = TRUE)
}

load_summary <- function(summary_rds = "output/results/mc_summary.rds",
                         summary_csv = "output/results/mc_summary.csv") {
  if (file.exists(summary_rds)) {
    return(readRDS(summary_rds))
  }
  if (file.exists(summary_csv)) {
    return(read.csv(summary_csv, stringsAsFactors = FALSE))
  }
  stop("Could not find summary file. Expected either:\n  - ", summary_rds, "\n  - ", summary_csv)
}

# pretty estimator labels
add_labels <- function(tab) {
  tab$rho_f   <- factor(tab$rho, levels = sort(unique(tab$rho)))
  tab$p_label <- paste0("p = ", tab$p_select)
  if ("p_hat_mean" %in% names(tab)) {
    tab$p_label <- paste0(tab$p_label, "\n(p̂≈", sprintf("%.2f", tab$p_hat_mean), ")")
  }
  tab$err_lab <- factor(tab$err_family, levels = c("normal", "t3", "t5", "logistic"))
  
  tab$est_lab <- factor(
    tab$estimator,
    levels = c("selected_ols", "zero_impute", "heckman_probit", "heckman_logit", "heckman_lpm"),
    labels = c("Selected OLS", "Zero-impute", "Heckman (probit)", "Heckman (logit)", "Heckman (LPM)")
  )
  tab
}

# -------------------------
# Main report generator
# -------------------------
make_tables_figs <- function(summary_rds = "output/results/mc_summary.rds") {
  make_dirs_if_needed()
  
  tab <- load_summary(summary_rds = summary_rds)
  
  # -------------------------
  # (A) Diagnostics tables
  # -------------------------
  diag_cols <- intersect(
    c("p_hat_mean", "p_hat_sd", "uv_corr_mean", "uv_corr_sd", "fail_rate", "n_fail", "n_used"),
    names(tab)
  )
  
  if (length(diag_cols) > 0) {
    diag_table <- tab %>%
      distinct(err_family, p_select, rho, !!!rlang::syms(diag_cols)) %>%
      arrange(err_family, p_select, rho)
  } else {
    diag_table <- tab %>%
      distinct(err_family, p_select, rho) %>%
      arrange(err_family, p_select, rho)
  }
  
  write.csv(diag_table,
            "output/results/tables/diagnostics_by_scenario.csv",
            row.names = FALSE)
  
  # Worst-case calibration deviation table (scenario-level)
  if ("p_hat_mean" %in% names(tab)) {
    calib_check <- tab %>%
      distinct(err_family, p_select, rho, p_hat_mean) %>%
      mutate(abs_dev = abs(p_hat_mean - p_select)) %>%
      arrange(desc(abs_dev))
    
    write.csv(calib_check,
              "output/results/tables/calibration_check_worst_cases.csv",
              row.names = FALSE)
  }
  
  # Failure rates by estimator and family
  if ("fail_rate" %in% names(tab)) {
    fail_by_est <- tab %>%
      distinct(estimator, err_family, p_select, rho, fail_rate, n_fail, n_used) %>%
      group_by(estimator, err_family) %>%
      summarise(
        max_fail_rate = max(fail_rate, na.rm = TRUE),
        total_fail    = sum(n_fail, na.rm = TRUE),
        min_n_used    = min(n_used, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(max_fail_rate), estimator, err_family)
    
    write.csv(fail_by_est,
              "output/results/tables/failure_summary.csv",
              row.names = FALSE)
  }
  
  # -------------------------
  # (B) Main results tables for beta1 (x1)
  # -------------------------
  b1 <- tab %>% filter(coef_name == "x1")
  
  # Main table (compact): average across rho (useful for RMSE/coverage overview)
  b1_avg_over_rho <- b1 %>%
    group_by(err_family, p_select, estimator) %>%
    summarise(
      bias      = mean(bias, na.rm = TRUE),
      rmse      = mean(rmse, na.rm = TRUE),
      sign_fail = mean(sign_fail, na.rm = TRUE),
      coverage  = mean(coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(err_family, p_select, estimator)
  
  write.csv(b1_avg_over_rho,
            "output/results/tables/main_beta1_avg_over_rho.csv",
            row.names = FALSE)
  
  # Appendix table: full grid by rho
  b1_full <- b1 %>%
    select(err_family, p_select, rho, estimator, bias, rmse, sign_fail, coverage) %>%
    arrange(err_family, p_select, rho, estimator)
  
  write.csv(b1_full,
            "output/results/tables/app_beta1_full_by_rho.csv",
            row.names = FALSE)
  
  # NEW: Main table for main text: focus on extreme selection (rho = ±0.6)
  b1_pm <- b1 %>%
    filter(rho %in% c(-0.6, 0.6)) %>%
    select(err_family, p_select, rho, estimator, bias, rmse, sign_fail, coverage) %>%
    arrange(err_family, p_select, rho, estimator)
  
  write.csv(b1_pm,
            "output/results/tables/main_beta1_rho_pm06.csv",
            row.names = FALSE)
  
  # Optional: Wide version (handy for LaTeX tables)
  b1_pm_wide <- b1_pm %>%
    pivot_wider(
      names_from = rho,
      values_from = c(bias, rmse, sign_fail, coverage),
      names_glue = "{.value}_rho{rho}"
    )
  
  write.csv(b1_pm_wide,
            "output/results/tables/main_beta1_rho_pm06_wide.csv",
            row.names = FALSE)
  
  # x2: save only bias/RMSE (sign_fail and coverage in your summary are for beta1 only)
  x2 <- tab %>% filter(coef_name == "x2")
  if (nrow(x2) > 0) {
    x2_avg <- x2 %>%
      group_by(err_family, p_select, estimator) %>%
      summarise(
        bias = mean(bias, na.rm = TRUE),
        rmse = mean(rmse, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(err_family, p_select, estimator)
    
    write.csv(x2_avg,
              "output/results/tables/secondary_x2_avg_over_rho.csv",
              row.names = FALSE)
  }
  
  # -------------------------
  # (C) Plots (beta1)
  # -------------------------
  b1p <- add_labels(b1)
  
  # 1) sign_fail
  p_sign <- ggplot(b1p, aes(x = rho_f, y = sign_fail, colour = est_lab, group = est_lab)) +
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
  
  ggsave("output/figs/sign_fail_beta1_by_rho.png", p_sign, width = 9, height = 6, dpi = 300)
  
  # 2) bias
  p_bias <- ggplot(b1p, aes(x = rho_f, y = bias, colour = est_lab, group = est_lab)) +
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
  
  ggsave("output/figs/bias_beta1_by_rho.png", p_bias, width = 9, height = 6, dpi = 300)
  
  # 3) coverage
  p_cov <- ggplot(b1p, aes(x = rho_f, y = coverage, colour = est_lab, group = est_lab)) +
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
  
  ggsave("output/figs/coverage_beta1_by_rho.png", p_cov, width = 9, height = 6, dpi = 300)
  
  # 4) optional: failure rate plot
  p_fail <- NULL
  if ("fail_rate" %in% names(b1p)) {
    p_fail <- ggplot(b1p, aes(x = rho_f, y = fail_rate, colour = est_lab, group = est_lab)) +
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
    
    ggsave("output/figs/fail_rate_by_rho.png", p_fail, width = 9, height = 6, dpi = 300)
  }
  
  # 5) optional: calibration plot p_hat_mean vs p_select
  p_calib <- NULL
  if ("p_hat_mean" %in% names(tab)) {
    calib_points <- tab %>%
      distinct(err_family, p_select, rho, p_hat_mean) %>%
      mutate(
        rho_f   = factor(rho, levels = sort(unique(rho))),
        err_lab = factor(err_family, levels = c("normal", "t3", "t5", "logistic"))
      )
    
    p_calib <- ggplot(calib_points, aes(x = p_select, y = p_hat_mean)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_point() +
      facet_grid(err_lab ~ rho_f) +
      labs(
        title = "Calibration check: achieved p̂ vs target p",
        x = "Target selection rate p",
        y = "Achieved selection rate p̂"
      ) +
      theme_minimal()
    
    ggsave("output/figs/calibration_check.png", p_calib, width = 9, height = 6, dpi = 300)
  }
  
  invisible(list(
    table_main_beta1_avg_over_rho = b1_avg_over_rho,
    table_main_beta1_rho_pm06     = b1_pm,
    table_main_beta1_rho_pm06_wide= b1_pm_wide,
    table_app_beta1_full_by_rho   = b1_full,
    plots = list(sign = p_sign, bias = p_bias, coverage = p_cov, fail = p_fail, calib = p_calib)
  ))
}

# Run if executed as a script
if (sys.nframe() == 0L) {
  make_tables_figs()
}



