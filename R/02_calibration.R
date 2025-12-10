# Selection-rate calibration
# Purpose: for each scenario, choose selection intercept gamma0 so Pr(S=1) ~= target p
# Inputs: scenario row (incl. error family for v), gamma_slopes, distribution of Z
# Outputs: gamma0 value (and optionally achieved p_hat)
# Sanity checks: achieved p_hat close to target; monotonicity in gamma0


`%||%` <- function(a, b) if (!is.null(a)) a else b

Fv_from_scen <- function(scen) {
  # CDF of v depending on scenario's error family
  if (scen$err_family == "normal" || scen$err_family == "mixture") {
    return(pnorm)
  }
  if (scen$err_family == "t3") {
    return(function(x) pt(x, df = 3))
  }
  if (scen$err_family == "t5") {
    return(function(x) pt(x, df = 5))
  }
  stop("Unknown err_family: ", scen$err_family)
}

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

add_gamma0_to_grid <- function(grid, Z_noInt, gamma_slopes) {
  # Adds gamma0 column for each row of grid (scenario table)
  grid$gamma0 <- NA_real_
  
  for (i in seq_len(nrow(grid))) {
    scen <- grid[i, ]
    Fv <- Fv_from_scen(scen)
    
    grid$gamma0[i] <- calibrate_gamma0(
      Z_noInt = Z_noInt,
      gamma_slopes = gamma_slopes,
      target_p = scen$p_select,
      Fv = Fv
    )
  }
  
  grid
}
