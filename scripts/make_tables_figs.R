# Reporting script
# Goal:
#   - load saved results
#   - produce summary tables (bias/RMSE/sign-fail/coverage/SE diagnostics)
#   - produce and save plots to output/figs/ and tables to output/tables/

suppressPackageStartupMessages({
  library(ggplot2)
})

# -----------------------------
# Helpers
# -----------------------------

make_dirs_if_needed <- function() {
  if (!dir.exists("output")) dir.create("output")
  if (!dir.exists("output/results")) dir.create("output/results", recursive = TRUE)
  if (!dir.exists("output/figs")) dir.create("output/figs", recursive = TRUE)
  if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
}

assert_has_cols <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
}

read_summary <- function(path = "output/results/mc_summary.rds") {
  if (!file.exists(path)) {
    stop("Summary file not found: ", path,
         "\nRun run_full_mc.R first to generate output/results/mc_summary.rds")
  }
  tab <- readRDS(path)
  
  # Required columns for most plots
  req <- c("scenario_id","n","p_select","rho","err_family","estimator",
           "coef_name","bias","rmse","sign_fail","sd_emp")
  assert_has_cols(tab, req)
  
  # Make types consistent for plotting
  tab$rho <- as.numeric(tab$rho)
  tab$n <- as.integer(tab$n)
  tab$p_select <- as.numeric(tab$p_select)
  
  tab$err_family <- factor(tab$err_family, levels = unique(tab$err_family))
  tab$estimator <- factor(tab$estimator, levels = unique(tab$estimator))
  tab$coef_name <- as.character(tab$coef_name)
  
  tab
}

# write a small LaTeX table with booktabs (no extra packages required beyond booktabs)
write_tex_table <- function(df, file, caption = NULL, label = NULL) {
  lines <- c()
  if (!is.null(caption)) lines <- c(lines, sprintf("\\caption{%s}", caption))
  if (!is.null(label))   lines <- c(lines, sprintf("\\label{%s}", label))
  
  # Build tabular header
  cols <- names(df)
  align <- paste(rep("l", length(cols)), collapse = " ")
  lines <- c(lines, sprintf("\\begin{tabular}{%s}", align))
  lines <- c(lines, "\\toprule")
  lines <- c(lines, paste(cols, collapse = " & "), " \\\\")
  lines <- c(lines, "\\midrule")
  
  # Body
  for (i in seq_len(nrow(df))) {
    row <- vapply(df[i, , drop = FALSE], function(x) as.character(x), character(1))
    lines <- c(lines, paste(row, collapse = " & "), " \\\\")
  }
  
  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}")
  
  # Wrap in table env for convenience
  wrapped <- c("\\begin{table}[!ht]", "\\centering",
               lines,
               "\\end{table}")
  
  writeLines(wrapped, con = file)
}

# -----------------------------
# Plotting functions
# -----------------------------

plot_metric_vs_rho <- function(df, ycol, ylab, title,
                               filename_base,
                               ylim01 = FALSE,
                               hline = NULL) {
  
  p <- ggplot(df, aes(x = rho, y = .data[[ycol]], group = estimator)) +
    geom_line(aes(linetype = estimator)) +
    geom_point(aes(shape = estimator)) +
    facet_grid(err_family ~ p_select) +
    labs(x = expression(rho), y = ylab, title = title) +
    theme_bw()
  
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = 2)
  }
  if (ylim01) p <- p + coord_cartesian(ylim = c(0, 1))
  
  # Save both pdf (paper) and png (quick view)
  ggsave(file.path("output/figs", paste0(filename_base, ".pdf")),
         plot = p, width = 10, height = 6)
  ggsave(file.path("output/figs", paste0(filename_base, ".png")),
         plot = p, width = 10, height = 6, dpi = 300)
  
  p
}

plot_se_ratio <- function(df, filename_base = "se_ratio_x1") {
  assert_has_cols(df, c("mean_se","sd_emp"))
  
  df$se_ratio <- df$mean_se / df$sd_emp
  
  p <- ggplot(df, aes(x = rho, y = se_ratio, group = estimator)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_line(aes(linetype = estimator)) +
    geom_point(aes(shape = estimator)) +
    facet_grid(err_family ~ p_select) +
    labs(x = expression(rho), y = "mean(SE) / sd(beta_hat)",
         title = paste0("SE calibration diagnostic for ", unique(df$coef_name))) +
    theme_bw()
  
  ggsave(file.path("output/figs", paste0(filename_base, ".pdf")),
         plot = p, width = 10, height = 6)
  ggsave(file.path("output/figs", paste0(filename_base, ".png")),
         plot = p, width = 10, height = 6, dpi = 300)
  
  p
}

plot_signfail_heatmap <- function(df, filename_base = "signfail_heatmap_x1") {
  # Heatmap: err_family (rows) x rho (cols), faceted by p_select and estimator
  df$rho_f <- factor(df$rho, levels = sort(unique(df$rho)))
  df$p_select <- factor(df$p_select, levels = sort(unique(df$p_select)))
  
  p <- ggplot(df, aes(x = rho_f, y = err_family, fill = sign_fail)) +
    geom_tile() +
    facet_grid(p_select ~ estimator) +
    labs(x = expression(rho), y = "Error family", fill = "Sign-fail",
         title = paste0("Sign-failure probability for ", unique(df$coef_name))) +
    theme_bw()
  
  ggsave(file.path("output/figs", paste0(filename_base, ".pdf")),
         plot = p, width = 11, height = 5)
  ggsave(file.path("output/figs", paste0(filename_base, ".png")),
         plot = p, width = 11, height = 5, dpi = 300)
  
  p
}

# Optional: selection calibration diagnostic using sel_rate_diag + mc_summary keys
plot_selection_diagnostic <- function(summary_tab,
                                      sel_path = "output/results/sel_rate_diag.rds",
                                      filename_base = "selection_rate_diag") {
  if (!file.exists(sel_path)) {
    message("Selection diagnostic file not found: ", sel_path, " (skipping selection diagnostic plot)")
    return(invisible(NULL))
  }
  sel <- readRDS(sel_path)
  assert_has_cols(sel, c("scenario_id","sel_rate_mean"))
  
  # pull scenario-level info from summary (one row per scenario_id)
  scen_info <- unique(summary_tab[, c("scenario_id","p_select","rho","err_family","n")])
  diag <- merge(scen_info, sel, by = "scenario_id", all.x = TRUE)
  
  p <- ggplot(diag, aes(x = p_select, y = sel_rate_mean)) +
    geom_point() +
    facet_grid(err_family ~ rho) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    labs(x = "Target selection rate (p_select)",
         y = "Mean realized selection rate",
         title = "Selection-rate calibration diagnostic") +
    theme_bw()
  
  ggsave(file.path("output/figs", paste0(filename_base, ".pdf")),
         plot = p, width = 10, height = 6)
  ggsave(file.path("output/figs", paste0(filename_base, ".png")),
         plot = p, width = 10, height = 6, dpi = 300)
  
  p
}

# -----------------------------
# Tables
# -----------------------------

# 1) Takeaway table: by p_select and estimator, show max/mean sign-fail and rmse over scenarios
make_takeaway_table <- function(df, focus_coef = "x1") {
  d <- subset(df, coef_name == focus_coef)
  
  # aggregate over all scenarios within each (p_select, estimator)
  agg_fun <- function(x) c(max = max(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE))
  
  a1 <- aggregate(sign_fail ~ p_select + estimator, data = d, FUN = agg_fun)
  a2 <- aggregate(rmse      ~ p_select + estimator, data = d, FUN = agg_fun)
  
  out <- merge(a1, a2, by = c("p_select","estimator"), all = TRUE)
  
  out$max_signfail  <- out$sign_fail[, "max"]
  out$mean_signfail <- out$sign_fail[, "mean"]
  out$max_rmse      <- out$rmse[, "max"]
  out$mean_rmse     <- out$rmse[, "mean"]
  out$sign_fail <- NULL
  out$rmse <- NULL
  
  # pretty formatting
  out$p_select <- sprintf("%.2f", as.numeric(out$p_select))
  out$max_signfail  <- sprintf("%.3f", out$max_signfail)
  out$mean_signfail <- sprintf("%.3f", out$mean_signfail)
  out$max_rmse      <- sprintf("%.3f", out$max_rmse)
  out$mean_rmse     <- sprintf("%.3f", out$mean_rmse)
  
  # order rows nicely
  out <- out[order(out$p_select, out$estimator), ]
  out
}

# 2) Optional: family-level table (more detailed, good for appendix)
make_family_table <- function(df, focus_coef = "x1") {
  d <- subset(df, coef_name == focus_coef)
  agg_fun <- function(x) c(mean = mean(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
  
  a1 <- aggregate(sign_fail ~ p_select + err_family + estimator, data = d, FUN = agg_fun)
  a2 <- aggregate(rmse      ~ p_select + err_family + estimator, data = d, FUN = agg_fun)
  
  out <- merge(a1, a2, by = c("p_select","err_family","estimator"), all = TRUE)
  
  out$mean_signfail <- out$sign_fail[, "mean"]
  out$max_signfail  <- out$sign_fail[, "max"]
  out$mean_rmse     <- out$rmse[, "mean"]
  out$max_rmse      <- out$rmse[, "max"]
  
  out$sign_fail <- NULL
  out$rmse <- NULL
  
  out
}

# -----------------------------
# Main driver
# -----------------------------

make_report <- function(summary_path = "output/results/mc_summary.rds",
                        focus_coef = "x1",
                        include_heatmap = TRUE,
                        include_sel_diag = TRUE) {
  make_dirs_if_needed()
  
  tab <- read_summary(summary_path)
  
  # Focus coefficient
  tab_f <- subset(tab, coef_name == focus_coef)
  
  # Core plots (paper-ready)
  plot_metric_vs_rho(tab_f, "bias", "Bias",
                     title = paste0("Bias for ", focus_coef),
                     filename_base = paste0("bias_", focus_coef),
                     hline = 0)
  
  plot_metric_vs_rho(tab_f, "rmse", "RMSE",
                     title = paste0("RMSE for ", focus_coef),
                     filename_base = paste0("rmse_", focus_coef))
  
  if ("coverage" %in% names(tab_f)) {
    plot_metric_vs_rho(tab_f, "coverage", "Coverage (95\\% CI)",
                       title = paste0("Coverage for ", focus_coef),
                       filename_base = paste0("coverage_", focus_coef),
                       ylim01 = TRUE, hline = 0.95)
  }
  
  plot_metric_vs_rho(tab_f, "sign_fail", "Sign-fail probability",
                     title = paste0("Sign-fail probability for ", focus_coef),
                     filename_base = paste0("signfail_", focus_coef),
                     ylim01 = TRUE, hline = 0.5)
  
  if (all(c("mean_se","sd_emp") %in% names(tab_f))) {
    plot_se_ratio(tab_f, filename_base = paste0("se_ratio_", focus_coef))
  }
  
  if (isTRUE(include_heatmap)) {
    plot_signfail_heatmap(tab_f, filename_base = paste0("signfail_heatmap_", focus_coef))
  }
  
  if (isTRUE(include_sel_diag)) {
    plot_selection_diagnostic(tab, sel_path = "output/results/sel_rate_diag.rds",
                              filename_base = "selection_rate_diag")
  }
  
  # Tables
  take <- make_takeaway_table(tab, focus_coef = focus_coef)
  write.csv(take, file = file.path("output/tables", paste0("takeaways_", focus_coef, ".csv")),
            row.names = FALSE)
  escape_latex <- function(x) {
    x <- as.character(x)
    x <- gsub("\\\\", "\\\\textbackslash{}", x)  # backslash first
    x <- gsub("([%&#_$])", "\\\\\\1", x, perl = TRUE)
    x <- gsub("~", "\\\\textasciitilde{}", x)
    x <- gsub("\\^", "\\\\textasciicircum{}", x)
    x
  }
  
  write_tex_table <- function(df, file, caption = NULL, label = NULL,
                              wrap_table_env = TRUE) {
    
    # escape column names + values
    colnames(df) <- escape_latex(names(df))
    df[] <- lapply(df, escape_latex)
    
    cols <- names(df)
    align <- paste(rep("l", length(cols)), collapse = "")
    
    lines <- c()
    lines <- c(lines, sprintf("\\begin{tabular}{%s}", align))
    lines <- c(lines, "\\toprule")
    lines <- c(lines, paste(cols, collapse = " & "), " \\\\")
    lines <- c(lines, "\\midrule")
    
    for (i in seq_len(nrow(df))) {
      row <- vapply(df[i, , drop = FALSE], identity, character(1))
      lines <- c(lines, paste(row, collapse = " & "), " \\\\")
    }
    
    lines <- c(lines, "\\bottomrule")
    lines <- c(lines, "\\end{tabular}")
    
    if (wrap_table_env) {
      wrapped <- c("\\begin{table}[!ht]", "\\centering",
                   if (!is.null(caption)) sprintf("\\caption{%s}", caption) else NULL,
                   if (!is.null(label)) sprintf("\\label{%s}", label) else NULL,
                   lines,
                   "\\end{table}")
      writeLines(wrapped, con = file)
    } else {
      writeLines(lines, con = file)
    }
  }
  
  fam <- make_family_table(tab, focus_coef = focus_coef)
  write.csv(fam, file = file.path("output/tables", paste0("by_family_", focus_coef, ".csv")),
            row.names = FALSE)
  # (Optional) tex for appendix; comment out if you prefer
  write_tex_table(fam,
                  file = file.path("output/tables", paste0("by_family_", focus_coef, ".tex")),
                  caption = paste0("Performance by error family for ", focus_coef, "."),
                  label = paste0("tab:byfamily_", focus_coef))
  
  message("Reporting complete.")
  invisible(list(summary = tab, takeaways = take, by_family = fam))
}

# If you want this to run when sourcing:
# make_report("output/results/mc_summary.rds", focus_coef = "x1")



