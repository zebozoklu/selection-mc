# Selection-rate calibration
# Purpose: for each scenario, choose selection intercept gamma0 so Pr(S=1) ~= target p
# Inputs: scenario row (incl. error family for v), gamma_slopes, distribution of Z
# Outputs: gamma0 value (and optionally achieved p_hat)
# Sanity checks: achieved p_hat close to target; monotonicity in gamma0
