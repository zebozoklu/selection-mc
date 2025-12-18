# Selection-rate calibration
# Purpose: for each scenario, choose selection intercept gamma0 so Pr(S=1) ~= target p
# Inputs: scenario row (incl. sel_model), gamma_slopes, distribution of Z
# Outputs: gamma0 value (and optionally achieved p_hat)
# Sanity checks: achieved p_hat close to target; monotonicity in gamma0

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Link function for selection probability given scenario
# Returns g(eta) = P(S=1 | eta), where eta = Z_noInt %*% gamma_slopes + gamma0
link_from_scen <- function(scen) {
  sel_model <- as.character(scen$sel_model)
  
  if (sel_model == "probit") {
    # probit: P(S=1|Z) = Phi(eta)
    return(function(eta) pnorm(eta))
  }
  
  if (sel_model == "logit") {
    # logit: P(S=1|Z) = 1/(1 + exp(-eta))
    return(function(eta) 1 / (1 + exp(-eta)))
  }
  
  if (sel_model == "lpm") {
    # linear probability model: clamp(eta, 0, 1)
    return(function(eta) {
      p <- eta
      p <- pmin(pmax(p, 0), 1)   # clamp into [0,1]
      p
    })
  }
  
  stop("Unknown sel_model: ", sel_model)
}

# Given Z_noInt and gamma_slopes, choose gamma0 so average selection prob ≈ target_p
calibrate_gamma0 <- function(Z_noInt, gamma_slopes, target_p, link_fun,
                             bracket = c(-10, 10)) {
  
  eta_noInt <- as.numeric(Z_noInt %*% gamma_slopes)
  
  f <- function(g0) {
    eta <- g0 + eta_noInt
    mean(link_fun(eta)) - target_p
  }
  
  lo <- bracket[1]; hi <- bracket[2]
  # expand bracket until it straddles 0
  while (f(lo) > 0) lo <- lo - 5
  while (f(hi) < 0) hi <- hi + 5
  
  uniroot(f, lower = lo, upper = hi)$root
}

# Adds gamma0 column for each row of grid (scenario table)
add_gamma0_to_grid <- function(grid, Z_noInt, gamma_slopes) {
  grid$gamma0 <- NA_real_
  
  for (i in seq_len(nrow(grid))) {
    scen <- grid[i, ]
    link_fun <- link_from_scen(scen)
    
    grid$gamma0[i] <- calibrate_gamma0(
      Z_noInt      = Z_noInt,
      gamma_slopes = gamma_slopes,
      target_p     = scen$p_select,
      link_fun     = link_fun
    )
  }
  
  grid
}
