# R/07_plots.R
# Plots for Monte Carlo results
# Inputs:
#   output/results/mc_summary.rds
#   output/results/sel_rate_diag.rds  (optional but recommended)
# Outputs:
#   output/figs/*.pdf and *.png
#
# Assumes mc_summary includes columns (from 06_summaries.R update):
#   mean, bias, rmse, sd_emp, mcse_mean,
#   coverage, coverage_mcse,
#   sign_fail, sign_fail_mcse,
#   mean_se, beta_true, n_ok_beta, n_ok_cov, etc.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

make_dirs_if_needed <- function() {
  dir.create("output", recursive = TRUE, showWarnings = FALSE)
  dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
  dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
}

.save_both <- function(p, stem, w = 10, h = 6, dpi = 300) {
  ggsave(file.path("output/figs", paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path("output/figs", paste0(stem, ".png")), p, width = w, height = h, dpi = dpi)
  invisible(p)
}

.pretty_estimators <- function(df) {
  # Keep original keys but show nicer labels
  map <- c(
    selected_ols  = "Selected OLS",
    zero_impute   = "Zero-impute OLS",
    heckman_2step = "Heckman 2-step"
  )
  df$estimator_pretty <- ifelse(df$estimator %in% names(map),
                                unname(map[df$estimator]),
                                df$estimator)
  # stable ordering
  ord <- unique(df$estimator_pretty)
  df$estimator_pretty <- factor(df$estimator_pretty, levels = ord)
  df
}

.facet_specs <- function(df) {
  # If multiple n values exist, include n on rows with err_family
  if (length(unique(df$n)) > 1L) {
    list(f = facet_grid(rows = vars(err_family, n), cols = vars(p_select)))
  } else {
    list(f = facet_grid(rows = vars(err_family), cols = vars(p_select)))
  }
}

.plot_metric_vs_rho <- function(df, ycol, ylab, stem,
                                ylims01 = FALSE,
                                add_hline = NULL,
                                use_pretty_est = TRUE) {
  
  if (use_pretty_est && !("estimator_pretty" %in% names(df))) {
    df <- .pretty_estimators(df)
    est_aes <- "estimator_pretty"
  } else {
    est_aes <- "estimator"
  }
  
  facet_obj <- .facet_specs(df)$f
  
  p <- ggplot(df, aes(x = rho, y = .data[[ycol]], group = .data[[est_aes]])) +
    geom_line(aes(linetype = .data[[est_aes]])) +
    geom_point(aes(shape = .data[[est_aes]])) +
    facet_obj +
    labs(x = expression(rho), y = ylab) +
    theme_bw()
  
  if (!is.null(add_hline)) {
    p <- p + geom_hline(yintercept = add_hline, linetype = "dashed")
  }
  if (ylims01) p <- p + coord_cartesian(ylim = c(0, 1))
  
  .save_both(p, stem)
}

.plot_metric_with_mcse <- function(df, ycol, mcse_col, ylab, stem,
                                   ylims01 = FALSE,
                                   add_hline = NULL,
                                   use_pretty_est = TRUE) {
  
  # Requires columns present
  if (!(ycol %in% names(df)) || !(mcse_col %in% names(df))) {
    warning("Skipping ", stem, " because columns missing: ", ycol, " or ", mcse_col)
    return(invisible(NULL))
  }
  
  if (use_pretty_est && !("estimator_pretty" %in% names(df))) {
    df <- .pretty_estimators(df)
    est_aes <- "estimator_pretty"
  } else {
    est_aes <- "estimator"
  }
  
  facet_obj <- .facet_specs(df)$f
  
  # 95% MCSE bars around the plotted statistic
  df <- df %>%
    mutate(
      lo = .data[[ycol]] - 1.96 * .data[[mcse_col]],
      hi = .data[[ycol]] + 1.96 * .data[[mcse_col]]
    )
  
  p <- ggplot(df, aes(x = rho, y = .data[[ycol]], group = .data[[est_aes]])) +
    geom_line(aes(linetype = .data[[est_aes]])) +
    geom_point(aes(shape = .data[[est_aes]])) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.02, alpha = 0.8) +
    facet_obj +
    labs(x = expression(rho), y = ylab) +
    theme_bw()
  
  if (!is.null(add_hline)) {
    p <- p + geom_hline(yintercept = add_hline, linetype = "dashed")
  }
  if (ylims01) p <- p + coord_cartesian(ylim = c(0, 1))
  
  .save_both(p, stem)
}

make_plots <- function(summary_path = "output/results/mc_summary.rds",
                       sel_diag_path = "output/results/sel_rate_diag.rds",
                       focus_coef = "x1") {
  
  make_dirs_if_needed()
  
  if (!file.exists(summary_path)) {
    stop("Summary file not found: ", summary_path,
         "\nRun run_full_mc.R first.")
  }
  
  tab <- readRDS(summary_path)
  
  # Basic cleaning
  tab <- tab %>%
    mutate(
      rho = as.numeric(rho),
      p_select = factor(p_select, levels = sort(unique(p_select))),
      err_family = factor(err_family, levels = unique(err_family)),
      estimator = factor(estimator, levels = unique(estimator)),
      n = factor(n, levels = sort(unique(n)))
    )
  
  # Focus on one coefficient (default x1)
  tab_f <- tab %>% filter(coef_name == focus_coef)
  
  # ---- 1) Bias with MCSE bars (MCSE of mean = MCSE of bias)
  .plot_metric_with_mcse(
    df = tab_f,
    ycol = "bias",
    mcse_col = "mcse_mean",
    ylab = paste0("Bias (", focus_coef, ")"),
    stem = paste0("bias_", focus_coef, "_vs_rho"),
    ylims01 = FALSE
  )
  
  # ---- 2) RMSE (no MCSE in table)
  .plot_metric_vs_rho(
    df = tab_f,
    ycol = "rmse",
    ylab = paste0("RMSE (", focus_coef, ")"),
    stem = paste0("rmse_", focus_coef, "_vs_rho"),
    ylims01 = FALSE
  )
  
  # ---- 3) Coverage with MCSE bars + dashed 0.95
  if ("coverage" %in% names(tab_f)) {
    .plot_metric_with_mcse(
      df = tab_f,
      ycol = "coverage",
      mcse_col = "coverage_mcse",
      ylab = paste0("Coverage (95% CI) — ", focus_coef),
      stem = paste0("coverage_", focus_coef, "_vs_rho"),
      ylims01 = TRUE,
      add_hline = 0.95
    )
  }
  
  # ---- 4) Sign-fail with MCSE bars
  # sign_fail is stored only on the target coefficient row (non-NA)
  df_sf <- tab %>% filter(!is.na(sign_fail))
  
  if (nrow(df_sf) > 0) {
    df_sf <- df_sf %>%
      mutate(
        rho = as.numeric(rho),
        p_select = factor(p_select, levels = sort(unique(p_select))),
        err_family = factor(err_family, levels = unique(err_family)),
        estimator = factor(estimator, levels = unique(estimator)),
        n = factor(n, levels = sort(unique(n)))
      )
    
    target_name <- unique(df_sf$coef_name)
    target_name <- if (length(target_name) == 1L) target_name else "target"
    
    .plot_metric_with_mcse(
      df = df_sf,
      ycol = "sign_fail",
      mcse_col = "sign_fail_mcse",
      ylab = paste0("Sign-reversal probability — ", target_name),
      stem = paste0("signfail_", target_name, "_vs_rho"),
      ylims01 = TRUE
    )
  }
  
  # ---- 5) SE calibration: ratio mean_se / sd_emp (best “interpretation” diagnostic)
  if (all(c("mean_se", "sd_emp") %in% names(tab_f))) {
    df_ratio <- tab_f %>%
      mutate(se_ratio = mean_se / sd_emp)
    
    .plot_metric_vs_rho(
      df = df_ratio,
      ycol = "se_ratio",
      ylab = paste0("SE calibration: mean(SE) / SD(beta_hat) — ", focus_coef),
      stem = paste0("se_ratio_", focus_coef, "_vs_rho"),
      ylims01 = FALSE,
      add_hline = 1.0
    )
  }
  
  # ---- 6) Selection-rate diagnostic plot (mean realized vs target p_select)
  if (file.exists(sel_diag_path)) {
    sel_diag <- readRDS(sel_diag_path)
    
    scen_meta <- tab %>%
      distinct(scenario_id, n, p_select, rho, err_family)
    
    sel2 <- sel_diag %>%
      left_join(scen_meta, by = "scenario_id") %>%
      mutate(
        p_select_num = as.numeric(as.character(p_select)),
        rho = as.numeric(as.character(rho))
      )
    
    # Plot mean realized selection vs target selection rate
    if (length(unique(sel2$n)) > 1L) {
      facet_obj <- facet_grid(rows = vars(err_family, n), cols = vars(rho))
    } else {
      facet_obj <- facet_grid(rows = vars(err_family), cols = vars(rho))
    }
    
    p_sel <- ggplot(sel2, aes(x = p_select_num, y = sel_rate_mean)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      facet_obj +
      labs(x = "Target selection rate (p_select)",
           y = "Mean realized selection rate") +
      theme_bw()
    
    .save_both(p_sel, "sel_rate_diag")
  }
  
  message("Saved plots to output/figs/")
  invisible(TRUE)
}


