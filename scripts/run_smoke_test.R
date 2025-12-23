# Smoke test (no full Monte Carlo)
# Goal: run 1 scenario, 1 dataset:
#   - load config + design + DGP + estimators
#   - calibrate gamma0
#   - simulate one dataset
#   - run each estimator once
#   - print estimates + SEs + key diagnostics

source("R/00_config.R")
source("R/01_scenarios.R")
source("R/02_calibration.R")
source("R/03_dgp.R")
source("R/04_estimators.R")

cfg <- make_config()

# Calibration draw for Z slopes
set.seed(cfg$seed_calib)
Z_noInt <- matrix(
  rnorm(cfg$n_calib * length(cfg$gamma_slopes)),
  ncol = length(cfg$gamma_slopes)
)

grid  <- make_scenario_grid()
grid2 <- add_gamma0_to_grid(grid, Z_noInt, cfg$gamma_slopes)

# Example scenario (must match your new grid defaults)
# n=1000, p=0.6, rho=0.0, logistic errors (change as you like)
scen_example <- subset(
  grid2,
  n == 1000 & p_select == 0.6 &
    rho == 0.0 & err_family == "logistic"
)[1, ]

if (nrow(scen_example) == 0) stop("No scenario matched your subset() filter.")

set.seed(cfg$seed_mc)
dat <- simulate_one_dataset(scen_example, cfg)

cat("Scenario:\n")
print(scen_example)

cat("\nRealized selection rate (mean S):\n")
print(mean(dat$S))

cat("\nShare observed Y:\n")
print(mean(!is.na(dat$Y)))

cat("\nHead of X:\n")
print(head(dat$X))

cat("\nSummary of observed Y (with NAs for unselected):\n")
print(summary(dat$Y))

# --- run each estimator once ---
estimators <- list(
  selected_ols = est_selected_ols,
  zero_impute  = est_zero_impute_ols,
  heckman_2step = est_heckman_2step
)

cat("\n\nEstimator outputs (beta and robust SE):\n")
for (nm in names(estimators)) {
  cat("\n---", nm, "---\n")
  out <- estimators[[nm]](dat)
  print(out$beta)
  cat("SE:\n")
  print(out$se)
  if (!is.null(out$aux)) {
    cat("Aux:\n")
    print(out$aux)
  }
}
