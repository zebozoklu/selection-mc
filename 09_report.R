# R/09_report.R
# Friend-style reporting layer:
#   - per-scenario density overlays / histograms
#   - table slice export
#   - NEW: appendix multipanel pages (2x3 or 3x2) for all scenarios
#
# Inputs:
#   output/results/mc_summary.rds
#   output/draws/draws_<scenario_id>.rds   (saved by run_full_mc.R)
# Outputs:
#   output/figs/report/*.pdf and *.png
#   output/tables/report_table_slice_<coef>.csv
#
# Example:
#   source("R/09_report.R")
#   make_report(c("S004","S001","S013"), coef_name="x1", make_hists=TRUE)
#   make_appendix_density_pages(coef_name="x1", ncol=3, nrow=2)  # all scenarios

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ----------------------------
# Directories
# ----------------------------
make_dirs <- function() {
  dir.create("output/figs/report", recursive = TRUE, showWarnings = FALSE)
  dir.create("output/tables",     recursive = TRUE, showWarnings = FALSE)
  dir.create("output/results",    recursive = TRUE, showWarnings = FALSE)
  dir.create("output/draws",      recursive = TRUE, showWarnings = FALSE)
}

.save_both <- function(p, stem, w = 10, h = 6, dpi = 300) {
  ggsave(file.path("output/figs/report", paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path("output/figs/report", paste0(stem, ".png")), p, width = w, height = h, dpi = dpi)
  invisible(p)
}

# ----------------------------
# Labels
# ----------------------------
pretty_est <- function(x) {
  dplyr::recode(
    x,
    selected_ols  = "Selected OLS",
    zero_impute   = "Zero-impute OLS",
    heckman_2step = "Heckman 2-step",
    .default = x
  )
}

scenario_title <- function(meta, coef_name) {
  paste0(
    meta$scenario_id, ": ",
    meta$err_family,
    ", p=", sprintf("%.2f", meta$p_select),
    ", rho=", sprintf("%.2f", meta$rho)
  )
}

# ----------------------------
# Read draws for one scenario (long format)
# ----------------------------
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
  if (is.null(colnames(first_mat))) {
    stop("raw_beta matrices have NULL colnames. Fix R/05_mc_engine.R to set colnames.")
  }
  if (!(coef_name %in% colnames(first_mat))) {
    stop("coef_name '", coef_name, "' not found for ", scenario_id,
         ". Available: ", paste(colnames(first_mat), collapse = ", "))
  }
  
  long <- bind_rows(lapply(names(raw_beta), function(est) {
    v <- raw_beta[[est]][, coef_name]
    v <- v[is.finite(v)]
    data.frame(
      scenario_id = scenario_id,
      estimator_raw = est,
      estimator = pretty_est(est),
      beta_hat = v,
      stringsAsFactors = FALSE
    )
  }))
  
  meta <- data.frame(
    scenario_id = scenario_id,
    n = scen$n,
    p_select = scen$p_select,
    rho = scen$rho,
    err_family = as.character(scen$err_family),
    stringsAsFactors = FALSE
  )
  
  list(long = long, meta = meta)
}

# ----------------------------
# Single-scenario density overlay (friend-style)
# ----------------------------
plot_density_overlay <- function(scenario_id, coef_name = "x1",
                                 beta_true_val = NULL,
                                 out_stem = NULL,
                                 xlim_use = NULL,
                                 show_title = TRUE) {
  
  out <- read_draws_long(scenario_id, coef_name)
  long <- out$long
  meta <- out$meta
  
  ttl <- if (isTRUE(show_title)) scenario_title(meta, coef_name) else NULL
  
  # default x-limits
  if (is.null(xlim_use)) {
    rng <- range(long$beta_hat, finite = TRUE)
    pad <- 0.08 * diff(rng)
    xlim_use <- c(rng[1] - pad, rng[2] + pad)
  }
  
  p <- ggplot(long, aes(x = beta_hat, color = estimator)) +
    geom_density(linewidth = 0.9) +
    labs(
      title = ttl,
      x = "Estimate",
      y = "Density",
      color = NULL
    ) +
    coord_cartesian(xlim = xlim_use) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 10)
    )
  
  # ONE true line
  if (!is.null(beta_true_val) && is.finite(beta_true_val)) {
    p <- p + geom_vline(xintercept = beta_true_val,
                        linetype = "dashed",
                        linewidth = 0.6,
                        alpha = 0.7)
  }
  
  if (!is.null(out_stem)) {
    .save_both(p, out_stem, w = 10, h = 6)
  }
  p
}

# ----------------------------
# Single-scenario histograms (optional appendix)
# ----------------------------
plot_histograms <- function(scenario_id, coef_name = "x1",
                            beta_true_val = NULL,
                            out_stem = NULL) {
  
  out <- read_draws_long(scenario_id, coef_name)
  long <- out$long
  meta <- out$meta
  
  ttl <- paste0("Histograms - ", scenario_title(meta, coef_name))
  
  p <- ggplot(long, aes(x = beta_hat)) +
    geom_histogram(bins = 40) +
    facet_wrap(~ estimator, scales = "free_y") +
    labs(title = ttl, x = "Estimate", y = "Count") +
    theme_bw()
  
  if (!is.null(beta_true_val) && is.finite(beta_true_val)) {
    p <- p + geom_vline(xintercept = beta_true_val,
                        linetype = "dashed",
                        linewidth = 0.6,
                        alpha = 0.7)
  }
  
  if (is.null(out_stem)) out_stem <- paste0("hist_", scenario_id, "_", coef_name)
  .save_both(p, out_stem, w = 11, h = 6)
}

# ----------------------------
# Export: compact table slice
# ----------------------------
export_table_slice <- function(scenario_ids, coef_name = "x1",
                               summary_path = "output/results/mc_summary.rds") {
  
  if (!file.exists(summary_path)) stop("Missing summary: ", summary_path)
  
  tab <- readRDS(summary_path)
  
  slice <- tab %>%
    filter(.data$scenario_id %in% scenario_ids,
           .data$coef_name == coef_name) %>%
    mutate(estimator = pretty_est(as.character(.data$estimator))) %>%
    select(
      scenario_id, n, p_select, rho, err_family, estimator,
      beta_true, mean, bias, sd_emp, rmse, mean_se, coverage, sign_fail
    ) %>%
    arrange(err_family, n, p_select, rho, estimator)
  
  out_csv <- file.path("output/tables", paste0("report_table_slice_", coef_name, ".csv"))
  write.csv(slice, out_csv, row.names = FALSE)
  
  message("Wrote table slice: ", out_csv)
  invisible(slice)
}

# ----------------------------
# Main driver (selected scenarios)
# ----------------------------
make_report <- function(
    scenario_ids,
    coef_name = "x1",
    make_hists = FALSE,
    summary_path = "output/results/mc_summary.rds"
) {
  
  make_dirs()
  
  if (!file.exists(summary_path)) stop("Missing summary: ", summary_path)
  
  bt <- readRDS(summary_path) %>%
    filter(.data$coef_name == coef_name) %>%
    select(.data$scenario_id, .data$beta_true) %>%
    distinct()
  bt_map <- setNames(bt$beta_true, bt$scenario_id)
  
  for (sid in scenario_ids) {
    p <- plot_density_overlay(
      scenario_id = sid,
      coef_name = coef_name,
      beta_true_val = bt_map[[sid]],
      out_stem = paste0("density_", sid, "_", coef_name)
    )
    
    if (isTRUE(make_hists)) {
      plot_histograms(
        scenario_id = sid,
        coef_name = coef_name,
        beta_true_val = bt_map[[sid]],
        out_stem = paste0("hist_", sid, "_", coef_name)
      )
    }
  }
  
  export_table_slice(
    scenario_ids = scenario_ids,
    coef_name = coef_name,
    summary_path = summary_path
  )
  
  message("Saved report figures to output/figs/report/")
  invisible(TRUE)
}

# ----------------------------
# NEW: Appendix multipanel pages for many scenarios
# ----------------------------
make_appendix_density_pages <- function(
    scenario_ids = NULL,                 # NULL = infer from mc_summary
    coef_name = "x1",
    summary_path = "output/results/mc_summary.rds",
    ncol = 3, nrow = 2,                  # 3x2 or 2x3
    filter_rho = NULL,                   # e.g. "<0" or "== -0.6" (see below)
    page_prefix = "appendix_density",
    width_per_col = 4.2,                 # tweak sizing
    height_per_row = 3.2
) {
  
  make_dirs()
  
  if (!file.exists(summary_path)) stop("Missing summary: ", summary_path)
  tab <- readRDS(summary_path)
  
  # infer all scenario ids if not provided
  if (is.null(scenario_ids)) {
    scenario_ids <- tab %>% distinct(scenario_id) %>% arrange(scenario_id) %>% pull(scenario_id)
  }
  
  # optional filter on rho using a simple string rule
  # filter_rho examples:
  #   "<0"       keeps negative rho
  #   "== -0.6"  keeps exactly -0.6
  #   "!= 0"
  if (!is.null(filter_rho)) {
    meta <- tab %>% distinct(scenario_id, rho)
    rho_num <- as.numeric(as.character(meta$rho))
    keep <- rep(TRUE, nrow(meta))
    
    if (filter_rho == "<0") keep <- rho_num < 0
    if (filter_rho == ">0") keep <- rho_num > 0
    if (filter_rho == "!=0") keep <- rho_num != 0
    if (filter_rho == "== -0.6") keep <- abs(rho_num - (-0.6)) < 1e-12
    if (filter_rho == "== 0") keep <- abs(rho_num - 0) < 1e-12
    if (filter_rho == "== 0.6") keep <- abs(rho_num - 0.6) < 1e-12
    
    scenario_ids <- meta$scenario_id[keep]
    scenario_ids <- sort(unique(scenario_ids))
  }
  
  # map beta_true
  bt <- tab %>%
    filter(.data$coef_name == coef_name) %>%
    select(.data$scenario_id, .data$beta_true) %>%
    distinct()
  bt_map <- setNames(bt$beta_true, bt$scenario_id)
  
  per_page <- ncol * nrow
  n_pages <- ceiling(length(scenario_ids) / per_page)
  
  # We use patchwork if available; otherwise fall back to gridExtra
  have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
  
  for (pg in seq_len(n_pages)) {
    idx <- ((pg - 1) * per_page + 1):min(pg * per_page, length(scenario_ids))
    sids <- scenario_ids[idx]
    
    # build each panel, and share x-limits within the page
    # First: load all draws on the page to compute pooled x-range
    all_vals <- c()
    panels <- list()
    
    # collect values for range
    for (sid in sids) {
      out <- read_draws_long(sid, coef_name)
      all_vals <- c(all_vals, out$long$beta_hat)
    }
    rng <- range(all_vals, finite = TRUE)
    pad <- 0.06 * diff(rng)
    xlim_page <- c(rng[1] - pad, rng[2] + pad)
    
    # now create panels
    for (sid in sids) {
      panels[[sid]] <- plot_density_overlay(
        scenario_id = sid,
        coef_name = coef_name,
        beta_true_val = bt_map[[sid]],
        out_stem = NULL,
        xlim_use = xlim_page,
        show_title = TRUE
      ) +
        theme(
          legend.position = "none",
          axis.title = element_blank()
        )
    }
    
    # collect legend from first panel (keep bottom legend once per page)
    legend_plot <- plot_density_overlay(
      scenario_id = sids[1],
      coef_name = coef_name,
      beta_true_val = bt_map[[sids[1]]],
      out_stem = NULL,
      xlim_use = xlim_page,
      show_title = FALSE
    ) + theme(legend.position = "bottom")
    
    # extract legend if patchwork is available; otherwise just leave legends off
    if (have_patchwork) {
      library(patchwork)
      
      # Turn panels into patchwork grid
      grid_plot <- wrap_plots(panels, ncol = ncol, nrow = nrow)
      
      # Add a shared title and legend row
      page_title <- paste0("Appendix: Sampling distributions (", coef_name, ") — page ", pg, "/", n_pages)
      
      final <- grid_plot +
        plot_annotation(title = page_title) &
        theme(plot.title = element_text(size = 12))
      
      # add legend as a separate row below
      final2 <- final / legend_plot + plot_layout(heights = c(1, 0.12))
      
      stem <- sprintf("%s_page_%02d_%s", page_prefix, pg, coef_name)
      .save_both(final2, stem,
                 w = width_per_col * ncol,
                 h = height_per_row * nrow + 1.0)
    } else {
      # fallback: save each panel individually if patchwork isn't installed
      warning("patchwork not installed; saving individual panels instead of multi-panel pages.")
      for (sid in sids) {
        .save_both(panels[[sid]], paste0("density_", sid, "_", coef_name),
                   w = 6, h = 4)
      }
    }
  }
  
  message("Saved appendix pages to output/figs/report/")
  invisible(TRUE)
}
