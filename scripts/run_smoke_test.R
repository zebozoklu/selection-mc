# Smoke test (no full Monte Carlo)
# Goal:
#   - load config + design + DGP + estimators
#   - calibrate gamma0 using DGP-consistent X distribution
#   - simulate one dataset for a single scenario
#   - run each estimator once
#   - print basic diagnostics (selection rate, estimates, etc.)

source("R/00_config.R")
source("R/01_scenarios.R")
source("R/02_calibration.R")
source("R/03_dgp.R")
source("R/04_estimators.R")

cfg <- make_config()

# --- Calibration covariates should match DGP: X ~ N(0, SigmaX) ---
set.seed(cfg$seed_calib)

SigmaX <- make_SigmaX(ar = cfg$SigmaX_ar)
p <- length(cfg$beta)

Lx <- chol(SigmaX)
Zx <- matrix(rnorm(cfg$n_calib * p), nrow = cfg$n_calib, ncol = p)
X_calib <- Zx %*% Lx  # this is Z_noInt in calibration (same dimension as gamma_slopes)

grid  <- make_scenario_grid()
grid2 <- add_gamma0_to_grid(grid, Z_noInt = X_calib, gamma_slopes = cfg$gamma_slopes)

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

cat("\nRealized selection rate (mean S) and DGP diagnostics:\n")
cat("mean(S)    =", mean(dat$S), "\n")
if (!is.null(dat$p_hat))   cat("p_hat      =", dat$p_hat, "\n")
if (!is.null(dat$uv_corr)) cat("corr(u,v)  =", dat$uv_corr, "\n")

cat("\nHead of X:\n")
print(head(dat$X))

cat("\nSummary of observed Y (with NAs for unselected):\n")
print(summary(dat$Y))

# --- Run estimators once on this dataset ---
cat("\nEstimator outputs (one dataset):\n")

b_sel   <- est_selected_ols(dat, cfg)
b_zero  <- est_zero_impute_ols(dat, cfg)
b_h_pro <- est_heckman_probit(dat, cfg)
b_h_log <- est_heckman_logit(dat, cfg)
b_h_lpm <- est_heckman_lpm(dat, cfg)

cat("\nSelected-sample OLS:\n")
print(b_sel)

cat("\nZero-imputation OLS:\n")
print(b_zero)

cat("\nHeckman 2-step (probit first stage):\n")
print(b_h_pro)

cat("\nHeckman 2-step (logit first stage, prob->probit IMR):\n")
print(b_h_log)

cat("\nHeckman 2-step (LPM first stage, prob->probit IMR):\n")
print(b_h_lpm)



