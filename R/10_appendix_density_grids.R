# R/10_appendix_density_grids.R
# Purpose: Appendix-ready 2x3 density grids per error family
# Output: paper/figs/fig_density_grid_<family>_<coef>.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

make_dirs <- function() {
  dir.create("output/figs/report", recursive = TRUE, showWarnings = FALSE)
  dir.create("output/draws",      recursive = TRUE, showWarnings = FALSE)
  dir.create("paper/figs",        recursive = TRUE, showWarnings = FALSE)
}

pretty_est <- function(x) {
  dplyr::recode(
    x,
    selected_ols  = "Selected OLS",
    zero_impute   = "Zero-impute OLS",
    heckman_2step = "Heckman 2-step",
    .default = x
  )
}

read_draws_long <- function(scenario_id, coef_name = "x1",
                            draws_dir = "output/draws") {
  fp <- file.path(draws_dir, paste0("draws_", scenario_id, ".rds"))
  if (!file.exists(fp)) stop("Draws file not found: ", fp)
  
  obj <- readRDS(fp)
  scen <- obj$scen
  raw_beta <- obj$raw_beta
  
  if (is.null(raw_beta) || !is.list(raw_beta) || length(raw_beta) == 0) {
    stop("raw_beta missing/empty in: ", fp)
  }
  
  first_mat <- raw_beta[[1]]
  if (is.null(colnames(first_mat))) stop("raw_beta has NULL colnames in: ", fp)
  if (!(coef_name %in% colnames(first_mat))) {
    stop("coef '", coef_name, "' not found in: ", fp,
         ". Available: ", paste(colnames(first_mat), collapse = ", "))
  }
  
  long <- bind_rows(lapply(names(raw_beta), function(est) {
    v <- raw_beta[[est]][, coef_name]
    v <- v[is.finite(v)]
    data.frame(
      scenario_id = scenario_id,
      estimator   = pretty_est(est),
      beta_hat    = v,
      stringsAsFactors = FALSE
    )
  }))
  
  long$p_select   <- as.numeric(scen$p_select)
  long$rho        <- as.numeric(scen$rho)
  long$err_family <- as.character(scen$err_family)
  long$n          <- as.integer(scen$n)
  
  long
}

plot_density_grid_family <- function(scenario_ids,
                                     family,
                                     coef_name = "x1",
                                     true_value = 0.05,
                                     draws_dir = "output/draws",
                                     out_png = NULL,
                                     width = 12,
                                     height = 8.5,
                                     dpi = 300) {
  
  make_dirs()
  
  df <- bind_rows(lapply(scenario_ids, read_draws_long, coef_name = coef_name, draws_dir = draws_dir))
  df <- df %>% filter(err_family == family)
  
  if (nrow(df) == 0) stop("No rows for family='", family, "'. Check scenario_ids or err_family labels.")
  
  # palette + legend order (same as your single-plot function)
  pal <- c(
    "Selected OLS"    = "#1f77b4",
    "Zero-impute OLS" = "#d62728",
    "Heckman 2-step"  = "#2ca02c"
  )
  df$estimator <- factor(df$estimator, levels = names(pal))
  
  # facet labels
  df$rho_f <- factor(df$rho, levels = c(-0.6, 0, 0.6))
  df$p_f   <- factor(df$p_select, levels = c(0.3, 0.6))
  
  # shared x-range within a family (keeps panels comparable)
  qs <- quantile(df$beta_hat, probs = c(0.01, 0.99), na.rm = TRUE)
  pad <- 0.10 * diff(qs)
  xlim_use <- c(qs[1] - pad, qs[2] + pad)
  
  p <- ggplot(df, aes(x = beta_hat, color = estimator)) +
    geom_density(linewidth = 0.9) +
    geom_vline(xintercept = true_value, linetype = "dashed", linewidth = 0.8, color = "black") +
    scale_color_manual(values = pal, drop = FALSE) +
    coord_cartesian(xlim = xlim_use) +
    facet_grid(p_f ~ rho_f,
               labeller = labeller(
                 p_f = function(x) paste0("p_select = ", x),
                 rho_f = function(x) paste0("rho = ", x)
               )) +
    labs(
      title = paste0("Sampling distributions for ", coef_name, " (err = ", family, ", n = ", unique(df$n), ")"),
      x = "Estimate",
      y = "Density",
      color = "Estimator"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 11)
    )
  
  if (is.null(out_png)) {
    out_png <- file.path("paper/figs", paste0("fig_density_grid_", family, "_", coef_name, ".png"))
  }
  ggsave(out_png, p, width = width, height = height, dpi = dpi)
  out_png
}

# Convenience wrapper: make all 5 grids
make_all_density_grids <- function(coef_name = "x1", true_value = 0.05) {
  # all scenario ids in your design: S001..S030
  all_ids <- sprintf("S%03d", 1:30)
  
  families <- c("normal", "logistic", "lpm", "t5", "t3")
  
  paths <- lapply(families, function(fam) {
    plot_density_grid_family(
      scenario_ids = all_ids,
      family = fam,
      coef_name = coef_name,
      true_value = true_value
    )
  })
  names(paths) <- families
  message("Saved density grids:\n", paste(unlist(paths), collapse = "\n"))
  invisible(paths)
}
