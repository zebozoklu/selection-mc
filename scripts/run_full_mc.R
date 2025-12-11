# Full Monte Carlo runner over all scenarios
# Goal:
#   - build scenario grid and calibrate gamma0
#   - run MC for each scenario
#   - collect tidy summaries into a table
#   - save results into output/results/

source("R/00_config.R")
source("R/01_scenarios.R")
source("R/02_calibration.R")
source("R/03_dgp.R")
source("R/04_estimators.R")
source("R/06_summaries.R")
source("R/05_mc_engine.R")

cfg <- make_config()

# 1) Build scenario grid and calibrate gamma0
set.seed(cfg$seed_calib)
Z_noInt <- matrix(
  rnorm(cfg$n_calib * length(cfg$gamma_slopes)),
  ncol = length(cfg$gamma_slopes)
)

grid  <- make_scenario_grid()
grid2 <- add_gamma0_to_grid(grid, Z_noInt, cfg$gamma_slopes)

cat("Number of scenarios:", nrow(grid2), "\n")

# 2) Define estimators
estimators <- list(
  selected_ols   = est_selected_ols,
  zero_impute    = est_zero_impute_ols,
  heckman_2step  = est_heckman_2step
)

# 3) Choose R for main runs
# For a first full run, it's safer to use R_debug.
# Later you can switch to cfg$R_full when everything is stable.
R_main <- cfg$R_debug   # e.g. 50; change to cfg$R_full for main experiment

# 4) Loop over scenarios
all_rows <- list()

for (i in seq_len(nrow(grid2))) {
  scen <- grid2[i, ]
  
  cat(sprintf("Running scenario %s (%d/%d): n=%d, p=%.2f, rho=%.2f, err=%s\n",
              scen$scenario_id, i, nrow(grid2),
              scen$n, scen$p_select, scen$rho,
              as.character(scen$err_family)))
  
  # use different seed per scenario to decorrelate simulations
  seed_i <- cfg$seed_mc + i
  
  mc_res <- run_mc_one_scenario(
    scen       = scen,
    cfg        = cfg,
    estimators = estimators,
    R          = R_main,
    seed       = seed_i
  )
  
  tab_i <- mc_result_to_table(mc_res)
  all_rows[[i]] <- tab_i
}

summary_table <- do.call(rbind, all_rows)

# 5) Save results
if (!dir.exists("output/results")) {
  dir.create("output/results", recursive = TRUE)
}

saveRDS(summary_table, file = "output/results/mc_summary.rds")

# optional: also save as CSV for quick viewing
write.csv(summary_table, file = "output/results/mc_summary.csv",
          row.names = FALSE)

cat("Saved summary table to output/results/mc_summary.rds and .csv\n")
