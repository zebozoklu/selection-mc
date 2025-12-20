# run_full_mc.R
# Full Monte Carlo runner over all scenarios
# Goal:
#   - build scenario grid and calibrate gamma0
#   - run MC for each scenario
#   - collect tidy summaries into a table
#   - (NEW) save replication-level draws for density plots
#   - (NEW) optionally create friend-style density overlay plots
#   - save results into output/results/

source("R/00_config.R")
source("R/01_scenarios.R")
source("R/02_calibration.R")
source("R/03_dgp.R")
source("R/04_estimators.R")
source("R/05_mc_engine.R")
source("R/06_summaries.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

cfg <- make_config()

# ----------------------------
# Options for draw saving + density plots
# ----------------------------
save_rep_draws     <- TRUE                 # saves output/draws/draws_<scenario_id>.rds
make_density_plots <- TRUE                 # creates output/figs/density/*.pdf/.png
density_coef_name  <- "x1"                 # coefficient to plot
density_scenarios  <- NULL                 # NULL = all; or c("S001","S010",...)

# ----------------------------
# Make output dirs early
# ----------------------------
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
dir.create("output/draws",   recursive = TRUE, showWarnings = FALSE)
dir.create("output/figs/density", recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helper: save draws for one scenario
# ----------------------------
save_draws_one_scenario <- function(scen, mc_res) {
  # mc_res is output of run_mc_one_scenario()
  # expected to contain raw_beta/raw_se lists keyed by estimator name
  saveRDS(
    list(
      scen     = scen,
      raw_beta = mc_res$raw_beta,
      raw_se   = mc_res$raw_se,
      sel_rate = mc_res$sel_rate
    ),
    file = file.path("output/draws", paste0("draws_", scen$scenario_id, ".rds"))
  )
}

# ----------------------------
# Helper: density plot for one scenario (overlay by estimator)
# ----------------------------
plot_density_one_scenario <- function(scen, mc_res, tab_i,
                                      coef_name = "x1",
                                      out_dir = "output/figs/density") {
  # tab_i is mc_result_to_table(mc_res) for this scenario; contains beta_true
  bt <- tab_i %>%
    filter(coef_name == !!coef_name) %>%
    pull(beta_true)
  
  beta_true_val <- if (length(bt) >= 1L && is.finite(bt[1])) bt[1] else NA_real_
  
  raw_beta <- mc_res$raw_beta
  if (is.null(raw_beta) || !is.list(raw_beta) || length(raw_beta) == 0) return(invisible(NULL))
  
  # Ensure coef exists in the matrices
  first_mat <- raw_beta[[1]]
  if (is.null(colnames(first_mat))) {
    # If your engine doesn't set colnames, we cannot safely select "x1"
    return(invisible(NULL))
  }
  if (!(coef_name %in% colnames(first_mat))) return(invisible(NULL))
  
  long <- bind_rows(lapply(names(raw_beta), function(est) {
    v <- raw_beta[[est]][, coef_name]
    v <- v[is.finite(v)]
    data.frame(estimator = est, beta_hat = v, stringsAsFactors = FALSE)
  }))
  
  if (nrow(long) == 0) return(invisible(NULL))
  
  ttl <- paste0(
    "Sampling distribution of ", coef_name, " — ",
    scen$scenario_id,
    " (n=", scen$n,
    ", p=", sprintf("%.2f", scen$p_select),
    ", rho=", sprintf("%.2f", scen$rho),
    ", err=", as.character(scen$err_family), ")"
  )
  
  p <- ggplot(long, aes(x = beta_hat, group = estimator)) +
    geom_density() +
    labs(x = paste0("Estimate: ", coef_name), y = "Density", title = ttl) +
    theme_bw()
  
  if (is.finite(beta_true_val)) {
    p <- p + geom_vline(xintercept = beta_true_val, linetype = "dashed")
  }
  
  base <- paste0("density_", scen$scenario_id, "_", coef_name)
  ggsave(file.path(out_dir, paste0(base, ".pdf")), p, width = 10, height = 6)
  ggsave(file.path(out_dir, paste0(base, ".png")), p, width = 10, height = 6, dpi = 200)
  
  invisible(p)
}

# ----------------------------
# 1) Build scenario grid and calibrate gamma0
# ----------------------------
set.seed(cfg$seed_calib)

SigmaX <- matrix(c(1, 0.3,
                   0.3, 1), 2, 2)
Lx <- chol(SigmaX)

Z_noInt <- matrix(rnorm(cfg$n_calib * 2), ncol = 2) %*% Lx

grid  <- make_scenario_grid()
grid2 <- add_gamma0_to_grid(grid, Z_noInt, cfg$gamma_slopes)

cat("Number of scenarios:", nrow(grid2), "\n")

# ----------------------------
# 2) Define estimators
# ----------------------------
estimators <- list(
  selected_ols   = est_selected_ols,
  zero_impute    = est_zero_impute_ols,
  heckman_2step  = est_heckman_2step
)

# ----------------------------
# 3) Choose R for main runs
# ----------------------------
R_main <- cfg$R_full

# ----------------------------
# 4) Loop over scenarios
# ----------------------------
all_rows <- vector("list", nrow(grid2))

sel_diag <- data.frame(
  scenario_id   = grid2$scenario_id,
  sel_rate_mean = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(grid2))) {
  scen <- grid2[i, , drop = FALSE]
  
  # Optional: run density only for selected scenarios
  do_density_this <- isTRUE(make_density_plots) &&
    (is.null(density_scenarios) || (scen$scenario_id %in% density_scenarios))
  
  cat(sprintf(
    "Running scenario %s (%d/%d): n=%d, p=%.2f, rho=%.2f, err=%s\n",
    scen$scenario_id, i, nrow(grid2),
    scen$n, scen$p_select, scen$rho,
    as.character(scen$err_family)
  ))
  
  # different seed per scenario
  seed_i <- cfg$seed_mc + i
  
  mc_res <- run_mc_one_scenario(
    scen       = scen,
    cfg        = cfg,
    estimators = estimators,
    R          = R_main,
    seed       = seed_i
  )
  
  # selection diagnostic (robust to NA)
  sel_diag$sel_rate_mean[i] <- mean(mc_res$sel_rate, na.rm = TRUE)
  
  # summary rows
  tab_i <- mc_result_to_table(mc_res)
  all_rows[[i]] <- tab_i
  
  # (NEW) save replication-level draws for density plots
  if (isTRUE(save_rep_draws)) {
    save_draws_one_scenario(scen, mc_res)
  }
  
  # (NEW) density overlay plot per scenario
  if (isTRUE(do_density_this)) {
    # wrap in try so one bad scenario doesn't kill the whole run
    try(
      plot_density_one_scenario(
        scen  = scen,
        mc_res = mc_res,
        tab_i = tab_i,
        coef_name = density_coef_name,
        out_dir = "output/figs/density"
      ),
      silent = TRUE
    )
  }
}

summary_table <- do.call(rbind, all_rows)

# ----------------------------
# 5) Save results
# ----------------------------
saveRDS(summary_table, file = "output/results/mc_summary.rds")
write.csv(summary_table, file = "output/results/mc_summary.csv", row.names = FALSE)

saveRDS(sel_diag, file = "output/results/sel_rate_diag.rds")
write.csv(sel_diag, file = "output/results/sel_rate_diag.csv", row.names = FALSE)

cat("Saved summary table to output/results/mc_summary.rds and .csv\n")
cat("Saved selection diagnostics to output/results/sel_rate_diag.rds and .csv\n")
if (isTRUE(save_rep_draws)) cat("Saved replication draws to output/draws/\n")
if (isTRUE(make_density_plots)) cat("Saved density plots to output/figs/density/\n")


