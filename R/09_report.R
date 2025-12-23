# R/09_report.R
# Purpose: simple, paper-ready density overlays (PNG only)
# Inputs:
#   output/draws/draws_<scenario_id>.rds
# Outputs:
#   output/figs/report/density_<scenario_id>_<coef>.png
#   (optional) copies to paper/figs/<custom_name>.png
#
# Example:
#   source("R/09_report.R")
#   make_density_pngs(
#     scenario_ids = c("S004","S001","S007"),
#     coef_name = "x1",
#     true_value = 0.05,
#     copy_to_paper = TRUE,
#     paper_name_map = c(
#       S004="fig_density_benchmark_x1.png",
#       S001="fig_density_normal_stress_x1.png",
#       S007="fig_density_t3_stress_x1.png"
#     )
#   )

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
      estimator = pretty_est(est),
      beta_hat = v,
      stringsAsFactors = FALSE
    )
  }))
  
  meta <- list(
    scenario_id = scenario_id,
    n = scen$n,
    p_select = scen$p_select,
    rho = scen$rho,
    err_family = as.character(scen$err_family)
  )
  
  list(long = long, meta = meta)
}

plot_density_png <- function(scenario_id, coef_name = "x1",
                             true_value = 0.05,
                             out_dir = "output/figs/report",
                             file_stem = NULL) {
  out <- read_draws_long(scenario_id, coef_name)
  long <- out$long
  meta <- out$meta
  
  # Clear, high-contrast palette (easy to tell apart)
  pal <- c(
    "Selected OLS"    = "#1f77b4",  # blue
    "Zero-impute OLS" = "#d62728",  # red
    "Heckman 2-step"  = "#2ca02c"   # green
  )
  
  # Stable estimator order in legend
  long$estimator <- factor(long$estimator, levels = names(pal))
  
  title_txt <- paste0(
    "Sampling distribution of ", coef_name, " - ", meta$scenario_id,
    " (n=", meta$n,
    ", p=", sprintf("%.2f", meta$p_select),
    ", rho=", sprintf("%.2f", meta$rho),
    ", err=", meta$err_family, ")"
  )
  
  # Use trimmed x-range so tails don't ruin readability
  qs <- quantile(long$beta_hat, probs = c(0.01, 0.99), na.rm = TRUE)
  pad <- 0.10 * diff(qs)
  xlim_use <- c(qs[1] - pad, qs[2] + pad)
  
  p <- ggplot(long, aes(x = beta_hat, color = estimator)) +
    geom_density(linewidth = 1.1) +
    geom_vline(xintercept = true_value,
               linetype = "dashed", linewidth = 0.9, color = "black") +
    scale_color_manual(values = pal, drop = FALSE) +
    coord_cartesian(xlim = xlim_use) +
    labs(
      title = title_txt,
      x = "Estimate",
      y = "Density",
      color = "Estimator"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank()
    )
  
  if (is.null(file_stem)) file_stem <- paste0("density_", scenario_id, "_", coef_name)
  out_png <- file.path(out_dir, paste0(file_stem, ".png"))
  ggsave(out_png, p, width = 11, height = 7, dpi = 300)
  
  invisible(out_png)
}

make_density_pngs <- function(scenario_ids,
                              coef_name = "x1",
                              true_value = 0.05,
                              copy_to_paper = FALSE,
                              paper_name_map = NULL,
                              out_dir = "output/figs/report",
                              paper_dir = "paper/figs") {
  
  make_dirs()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (sid in scenario_ids) {
    stem <- paste0("density_", sid, "_", coef_name)
    out_png <- plot_density_png(
      scenario_id = sid,
      coef_name = coef_name,
      true_value = true_value,
      out_dir = out_dir,
      file_stem = stem
    )
    
    if (isTRUE(copy_to_paper)) {
      if (!is.null(paper_name_map) && sid %in% names(paper_name_map)) {
        dst <- file.path(paper_dir, unname(paper_name_map[[sid]]))
      } else {
        dst <- file.path(paper_dir, basename(out_png))
      }
      file.copy(out_png, dst, overwrite = TRUE)
    }
  }
  
  message("Saved PNGs to: ", out_dir)
  if (isTRUE(copy_to_paper)) message("Copied selected PNGs to: ", paper_dir)
  invisible(TRUE)
}



