# DGP layer
# Purpose: simulate one dataset from the selection model
# Model:
#   Y* = X beta + u
#   S  = 1{ Z gamma + v > 0 }
#   Y observed only if S=1
# Inputs: scenario + parameters (beta, gamma, rho, error family + params)
# Outputs: list(Y, S, X, Z, y_star, s_star)
# Sanity checks: mean(S) near target; dimensions consistent; no NA explosions
