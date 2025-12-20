# R/09_report.R
# Friend-style reporting layer:
#   - per-scenario density overlays / histograms
#   - table slice export
#   - appendix multipanel pages (2x3 or 3x2) for all scenarios
#   - (NEW) optional copy-to-paper/figs with stable filenames
#
# Inputs:
#   output/results/mc_summary.rds
#   output/draws/draws_<scenario_id>.rds   (saved by run_full_mc.R)
# Outputs:
#   output/figs/report/*.pdf and *.png
#   output/tables/report_table_slice_<coef>.csv
#   (optional) paper/figs/*.pdf (copied for LaTeX)
#
# Example:
#   source("R/09_report.R")
#   make_report(c("S004","S001","S007"), coef_name="x1", make_hists=FALSE,
#               copy_to_paper=TRUE)
#   make_appendix_density_pages(coef_name="x1", ncol=3, nrow=2)

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

.save_both <- function(p, stem, out_dir = "output/figs/report",
                       w = 10, h = 6, dpi = 300) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = dpi)
  invisible(p)
}

.copy_pdf_to_paper <- function(stem, out_dir = "output/figs/report",
                               paper_figs_dir = "paper/figs",
                               paper_filename = NULL) {
  dir.create(paper_figs_dir, recursive = TRUE, showWarnings = FALSE)
  src <- file.path(out_dir, paste0(stem, ".pdf"))
  if (!file.exists(src)) return(invisible(FALSE))
  
  if (is.null(paper_filename)) paper_filename <- paste0(stem, ".pdf")
  dst <- file.path(paper_figs_dir, paper_filename)
  
  ok <- file.copy(src, dst, overwrite = TRUE)
  invisible(ok)
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

# Friend-like title (avoid Unicode em dash to prevent LaTeX bookmark warnings)
scenario_title <- function(meta, coef_name) {
  paste0(
    "Sampling distribution of ", coef_name,
    " - ", meta$scenario_id,
    " (n=", meta$n,
    ", p=", sprintf("%.2f", meta$p_select),
    ", rho=", sprintf("%.2f", meta$rho),
    ", err=", as.character(meta$err_family), ")"
  )
}

# Manual palette to mimic friend-style contrast
# (You explicitly asked to match style, so fixed colors are OK.)
.est_palette <- function() {
  c(
    "Selected OLS"    = "grey50",
    "Zero-impute OLS" = "red3",
    "Heckman 2-step"  = "black"
  )
}

.friend_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 6, face = "plain", hjust = 0.5),
      axis.title = element_text(size = base_size + 2),
      axis.text  = element_text(size = base_size),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.box = "horizontal",
      panel.grid.minor = element_blank()
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
  
  # stable estimator order in legend
  long$estimator <- factor(long$estimator,
                           levels = c("Selected OLS", "Zero-impute OLS", "Heckman 2-step"))
  
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
                                 show_title = TRUE,
                                 base_size = 14) {
  
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
  
  pal <- .est_palette()
  
  p <- ggplot(long, aes(x = beta_hat, color = estimator)) +
    geom_density(linewidth = 1.1, adjust = 1.0) +
    # True value line (single)
    { if (!is.null(beta_true_val) && is.finite(beta_true_val))
      geom_vline(xintercept = beta_true_val, linetype = "dashed", linewidth = 0.9)
      else NULL } +
    scale_color_manual(values = pal, drop = FALSE) +
    coord_cartesian(xlim = xlim_use) +
    labs(
      title = ttl,
      x = "Estimate",
      y = "Density"
    ) +
    .friend_theme(base_size = base_size) +
    # Make legend keys look more like the friend example (square outlines instead of line segments)
    guides(
      color = guide_legend(
        override.aes = list(linetype = 0, shape = 0, size = 6, linewidth = 1.5)
      )
    )
  
  if (!is.null(out_stem)) {
    .save_both(p, out_stem, w = 11, h = 7)
  }
  p
}

# ----------------------------
# Single-scenario histograms (optional)
# ----------------------------
plot_histograms <- function(scenario_id, coef_name = "x1",
                            beta_true_val = NULL,
                            out_stem = NULL,
                            base_size = 13) {
  
  out <- read_draws_long(scenario_id, coef_name)
  long <- out$long
  meta <- out$meta
  
  ttl <- paste0("Histograms - ", scenario_title(meta, coef_name))
  
  p <- ggplot(long, aes(x = beta_hat, fill = estimator)) +
    geom_histogram(bins = 40, color = "white") +
    facet_wrap(~ estimator, scales = "free_y") +
    { if (!is.null(beta_true_val) && is.finite(beta_true_val))
      geom_vline(xintercept = beta_true_val, linetype = "dashed", linewidth = 0.9)
      else NULL } +
    labs(title = ttl, x = "Estimate", y = "Count") +
    theme_minimal(base_size = base_size) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  if (is.null(out_stem)) out_stem <- paste0("hist_", scenario_id, "_", coef_name)
  .save_both(p, out_stem, w = 11, h = 7)
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
    summary_path = "output/results/mc_summary.rds",
    copy_to_paper = FALSE,
    paper_figs_dir = "paper/figs",
    paper_name_map = NULL  # named character vector: c(S004="fig_density_benchmark_x1.pdf", ...)
) {
  
  make_dirs()
  if (!file.exists(summary_path)) stop("Missing summary: ", summary_path)
  
  bt <- readRDS(summary_path) %>%
    filter(.data$coef_name == coef_name) %>%
    select(.data$scenario_id, .data$beta_true) %>%
    distinct()
  bt_map <- setNames(bt$beta_true, bt$scenario_id)
  
  for (sid in scenario_ids) {
    stem <- paste0("density_", sid, "_", coef_name)
    
    plot_density_overlay(
      scenario_id = sid,
      coef_name = coef_name,
      beta_true_val = bt_map[[sid]],
      out_stem = stem,
      show_title = TRUE
    )
    
    if (isTRUE(make_hists)) {
      plot_histograms(
        scenario_id = sid,
        coef_name = coef_name,
        beta_true_val = bt_map[[sid]],
        out_stem = paste0("hist_", sid, "_", coef_name)
      )
    }
    
    if (isTRUE(copy_to_paper)) {
      paper_fn <- NULL
      if (!is.null(paper_name_map) && sid %in% names(paper_name_map)) {
        paper_fn <- unname(paper_name_map[[sid]])
      }
      .copy_pdf_to_paper(stem, paper_figs_dir = paper_figs_dir, paper_filename = paper_fn)
    }
  }
  
  export_table_slice(
    scenario_ids = scenario_ids,
    coef_name = coef_name,
    summary_path = summary_path
  )
  
  message("Saved report figures to output/figs/report/")
  if (isTRUE(copy_to_paper)) message("Copied selected PDFs to ", paper_figs_dir)
  invisible(TRUE)
}

# ----------------------------
# Appendix multipanel pages (many scenarios)
# ----------------------------
make_appendix_density_pages <- function(
    scenario_ids = NULL,                 # NULL = infer from mc_summary
    coef_name = "x1",
    summary_path = "output/results/mc_summary.rds",
    ncol = 3, nrow = 2,
    filter_rho = NULL,                   # "<0", "== -0.6", "!=0", "== 0", "== 0.6"
    page_prefix = "appendix_density",
    width_per_col = 4.4,
    height_per_row = 3.4
) {
  
  make_dirs()
  if (!file.exists(summary_path)) stop("Missing summary: ", summary_path)
  tab <- readRDS(summary_path)
  
  if (is.null(scenario_ids)) {
    scenario_ids <- tab %>% distinct(scenario_id) %>% arrange(scenario_id) %>% pull(scenario_id)
  }
  
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
    
    scenario_ids <- sort(unique(meta$scenario_id[keep]))
  }
  
  bt <- tab %>%
    filter(.data$coef_name == coef_name) %>%
    select(.data$scenario_id, .data$beta_true) %>%
    distinct()
  bt_map <- setNames(bt$beta_true, bt$scenario_id)
  
  per_page <- ncol * nrow
  n_pages <- ceiling(length(scenario_ids) / per_page)
  
  have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
  
  for (pg in seq_len(n_pages)) {
    idx <- ((pg - 1) * per_page + 1):min(pg * per_page, length(scenario_ids))
    sids <- scenario_ids[idx]
    
    # shared x-range within page
    all_vals <- c()
    for (sid in sids) {
      out <- read_draws_long(sid, coef_name)
      all_vals <- c(all_vals, out$long$beta_hat)
    }
    rng <- range(all_vals, finite = TRUE)
    pad <- 0.06 * diff(rng)
    xlim_page <- c(rng[1] - pad, rng[2] + pad)
    
    panels <- list()
    for (sid in sids) {
      panels[[sid]] <- plot_density_overlay(
        scenario_id = sid,
        coef_name = coef_name,
        beta_true_val = bt_map[[sid]],
        out_stem = NULL,
        xlim_use = xlim_page,
        show_title = TRUE,
        base_size = 10
      ) +
        theme(
          legend.position = "none",
          axis.title = element_blank(),
          plot.title = element_text(size = 11)
        )
    }
    
    if (have_patchwork) {
      library(patchwork)
      grid_plot <- wrap_plots(panels, ncol = ncol, nrow = nrow)
      page_title <- paste0("Appendix: Sampling distributions (", coef_name, ") - page ", pg, "/", n_pages)
      
      final <- grid_plot +
        plot_annotation(title = page_title) &
        theme(plot.title = element_text(size = 13, hjust = 0.5))
      
      stem <- sprintf("%s_page_%02d_%s", page_prefix, pg, coef_name)
      .save_both(final, stem,
                 w = width_per_col * ncol,
                 h = height_per_row * nrow)
    } else {
      warning("patchwork not installed; saving individual panels instead of multi-panel pages.")
      for (sid in sids) {
        .save_both(panels[[sid]], paste0("density_", sid, "_", coef_name),
                   w = 7, h = 5)
      }
    }
  }
  
  message("Saved appendix pages to output/figs/report/")
  invisible(TRUE)
}


