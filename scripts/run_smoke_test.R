# Smoke test (no full Monte Carlo)
# Goal: run 1 scenario, 1 replication:
#   - simulate one dataset
#   - run each estimator once
#   - print estimates + key diagnostics (selection rate, etc.)

# Smoke test: one scenario, one dataset
# Goal:
#   - load config + design + DGP code
#   - calibrate gamma0
#   - simulate one dataset
#   - print basic diagnostics (selection rate, head(X), summary(Y))

source("R/00_config.R")
source("R/01_scenarios.R")
source("R/02_calibration.R")
source("R/03_dgp.R")

cfg <- make_config()

set.seed(cfg$seed_calib)
Z_noInt <- matrix(
  rnorm(cfg$n_calib * length(cfg$gamma_slopes)),
  ncol = length(cfg$gamma_slopes)
)

grid  <- make_scenario_grid()
grid2 <- add_gamma0_to_grid(grid, Z_noInt, cfg$gamma_slopes)

# Example scenario: n=2000, p=0.6, rho=0.3, normal errors
scen_example <- subset(
  grid2,
  n == 2000 & p_select == 0.6 &
    rho == 0.3 & err_family == "normal"
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
