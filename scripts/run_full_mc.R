# Full Monte Carlo runner over all scenarios
# Goal:
#   - build scenario grid and calibrate gamma0 (DGP-consistent)
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

# 1) Build scenario grid and calibrate gamma0 using DGP-consistent X distribution
set.seed(cfg$seed_calib)

# X_calib ~ N(0, SigmaX) with the SAME SigmaX as in the DGP
SigmaX <- make_SigmaX(ar = cfg$SigmaX_ar)
p <- length(cfg$beta)

Lx <- chol(SigmaX)
Zx <- matrix(rnorm(cfg$n_calib * p), nrow = cfg$n_calib, ncol = p)
X_calib <- Zx %*% Lx   # this is Z_noInt for calibration

grid  <- make_scenario_grid()
grid2 <- add_gamma0_to_grid(grid, Z_noInt = X_calib, gamma_slopes = cfg$gamma_slopes)

cat("Number of scenarios:", nrow(grid2), "\n")

# 2) Define estimators
estimators <- list(
  selected_ols    = est_selected_ols,
  zero_impute     = est_zero_impute_ols,
  heckman_probit  = est_heckman_probit,
  heckman_logit   = est_heckman_logit,
  heckman_lpm     = est_heckman_lpm
)

# 3) Choose R for main runs
R_main <- cfg$R_full # change to cfg$R_full for main experiment

# 4) Loop over scenarios
all_rows <- vector("list", nrow(grid2))

for (i in seq_len(nrow(grid2))) {
  scen <- grid2[i, ]
  
  cat(sprintf(
    "Running scenario %s (%d/%d): n=%d, p=%.2f, rho=%.2f, err=%s\n",
    scen$scenario_id, i, nrow(grid2),
    scen$n, scen$p_select, scen$rho, as.character(scen$err_family)
  ))
  
  # use different seed per scenario to decorrelate simulations
  seed_i <- cfg$seed_mc + i
  
  mc_res <- run_mc_one_scenario(
    scen       = scen,
    cfg        = cfg,
    estimators = estimators,
    R          = R_main,
    seed       = seed_i
  )
  
  all_rows[[i]] <- mc_result_to_table(mc_res)
}

summary_table <- do.call(rbind, all_rows)

# 5) Save results
if (!dir.exists("output/results")) {
  dir.create("output/results", recursive = TRUE)
}

saveRDS(summary_table, file = "output/results/mc_summary.rds")

# optional: also save as CSV for quick viewing
write.csv(
  summary_table,
  file      = "output/results/mc_summary.csv",
  row.names = FALSE
)

cat("Saved summary table to output/results/mc_summary.rds and .csv\n")



