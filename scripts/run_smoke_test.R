# Smoke test (no full Monte Carlo)
# Goal:
#   - load config + design + DGP + estimators
#   - calibrate gamma0
#   - simulate one dataset for a single scenario
#   - run each estimator once
#   - print basic diagnostics (selection rate, estimates, etc.)

source("R/00_config.R")
source("R/01_scenarios.R")
source("R/02_calibration.R")
source("R/03_dgp.R")
source("R/04_estimators.R")

cfg <- make_config()

set.seed(cfg$seed_calib)
Z_noInt <- matrix(
  rnorm(cfg$n_calib * length(cfg$gamma_slopes)),
  ncol = length(cfg$gamma_slopes)
)

grid  <- make_scenario_grid()
grid2 <- add_gamma0_to_grid(grid, Z_noInt, cfg$gamma_slopes)

# Example scenario: n=500, p=0.3, rho=0.6, normal errors
scen_example <- subset(
  grid2,
  n == 500 &
    p_select == 0.3 &
    rho == 0.6 &
    err_family == "normal"
)[1, ]

set.seed(cfg$seed_mc)
dat <- simulate_one_dataset(scen_example, cfg)

cat("Scenario:\n")
print(scen_example)

cat("\nRealized selection rate (mean S):\n")
print(mean(dat$S))

cat("\nHead of X:\n")
print(head(dat$X))

cat("\nSummary of observed Y (with NAs for unselected):\n")
print(summary(dat$Y))

# --- Run estimators once on this dataset ---
cat("\nEstimator outputs (one dataset):\n")

b_sel  <- est_selected_ols(dat)
b_zero <- est_zero_impute_ols(dat)
b_heck <- est_heckman_2step(dat)

cat("\nSelected-sample OLS:\n")
print(b_sel)

cat("\nZero-imputation OLS:\n")
print(b_zero)

cat("\nHeckman 2-step (probit first stage):\n")
print(b_heck)



