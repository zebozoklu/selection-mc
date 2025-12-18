# Selection-rate calibration
# Purpose: for each scenario, choose selection intercept gamma0 so Pr(S=1) ~= target p
# Inputs: scenario row (incl. err_family), gamma_slopes, distribution of Z
# Outputs: gamma0 value (and optionally achieved p_hat)
# Sanity checks: achieved p_hat close to target; monotonicity in gamma0

`%||%` <- function(a, b) if (!is.null(a)) a else b

# CDF of v depending on scenario's error family
Fv_from_scen <- function(scen) {
  efam <- as.character(scen$err_family)
  
  if (efam == "normal") {
    return(pnorm)
  }
  if (efam == "t3") {
    return(function(x) pt(x, df = 3))
  }
  if (efam == "t5") {
    return(function(x) pt(x, df = 5))
  }
  if (efam == "logistic") {
    s <- sqrt(3) / pi
    return(function(x) plogis(x / s))
  }
  
  stop("Unknown err_family in Fv_from_scen: ", efam)
}

# Given Z_noInt and gamma_slopes, choose gamma0 so average selection prob ≈ target_p
calibrate_gamma0 <- function(Z_noInt, gamma_slopes, target_p, Fv = pnorm,
                             bracket = c(-10, 10)) {
  
  eta <- as.numeric(Z_noInt %*% gamma_slopes)
  
  f <- function(g0) mean(Fv(g0 + eta)) - target_p
  
  lo <- bracket[1]; hi <- bracket[2]
  # expand bracket until it straddles 0
  while (f(lo) > 0) lo <- lo - 5
  while (f(hi) < 0) hi <- hi + 5
  
  uniroot(f, lower = lo, upper = hi)$root
}

# Adds gamma0 column for each row of grid (scenario table)
add_gamma0_to_grid <- function(grid, Z_noInt, gamma_slopes) {
  grid$gamma0    <- NA_real_
  grid$p_calib   <- NA_real_   # NEW: achieved mean selection prob in calibration sample
  
  for (i in seq_len(nrow(grid))) {
    scen <- grid[i, ]
    Fv   <- Fv_from_scen(scen)
    
    g0 <- calibrate_gamma0(
      Z_noInt      = Z_noInt,
      gamma_slopes = gamma_slopes,
      target_p     = scen$p_select,
      Fv           = Fv
    )
    
    eta <- as.numeric(Z_noInt %*% gamma_slopes)
    p_calib <- mean(Fv(g0 + eta))
    
    grid$gamma0[i]  <- g0
    grid$p_calib[i] <- p_calib
  }
  
  grid
}





