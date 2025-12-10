# Estimator layer
# Purpose: implement each estimator mapping (Y,S,X,Z) -> beta_hat
# Estimators:
#   - selected-sample OLS
#   - zero-imputation OLS
#   - Heckman two-step (probit + IMR)
#   - optional series/spline correction
# Inputs: simulated dataset list
# Outputs: standardized beta_hat vector (intercept + slopes), plus optional extras
# Sanity checks: returns numeric vector with correct length; no silent NA
