# DGP layer
# Purpose: simulate one dataset from the selection model
# Model:
#   Y* = X beta + u
#   S  = 1{ Z gamma + v > 0 }
#   Y observed only if S=1
# Inputs: scenario + parameters (beta, gamma, rho, error family + params)
# Outputs: list(Y, S, X, Z, y_star, s_star)
# Sanity checks: mean(S) near target; dimensions consistent; no NA explosions


# --- helper: build SigmaX if not supplied (AR(0.3) for now) ---
make_SigmaX <- function(ar = 0.3) {
  matrix(c(1, ar,
           ar, 1), nrow = 2, ncol = 2)
}

# --- helper: draw bivariate normal errors (u, v) ---
draw_uv_normal <- function(n, rho) {
  Sigma <- matrix(c(1, rho,
                    rho, 1), 2, 2)
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * 2), n, 2)
  E <- Z %*% L
  list(u = E[,1], v = E[,2])
}

# --- helper: draw bivariate t errors via scaled normal ---
draw_uv_t <- function(n, rho, df) {
  base <- draw_uv_normal(n, rho)
  s <- sqrt(df / rchisq(n, df))
  list(
    u = base$u * s,
    v = base$v * s
  )
}

# --- helper: mixture case (currently: mixture in u only, v normal/t) ---
draw_uv_mixture <- function(n, rho, base_family = "normal", df = NULL,
                            mix_pi = 0.10, mix_shift_u = 2.0) {
  # base errors
  if (base_family == "normal") {
    base <- draw_uv_normal(n, rho)
  } else if (base_family == "t") {
    if (is.null(df)) stop("df must be provided for t mixture.")
    base <- draw_uv_t(n, rho, df = df)
  } else {
    stop("Unknown base_family in mixture: ", base_family)
  }
  
  jump <- rbinom(n, 1, mix_pi)
  u <- base$u + jump * mix_shift_u
  v <- base$v  # kept as base (normal/t) for now
  
  list(u = u, v = v)
}

simulate_one_dataset <- function(scen, cfg, SigmaX = NULL) {
  # scen: one-row data.frame with n, rho, err_family, df, mix_pi, mix_shift_u, gamma0
  # cfg: config list from make_config()
  
  # --- unpack basic stuff ---
  n    <- scen$n
  rho  <- scen$rho
  beta <- cfg$beta
  gamma_slopes <- cfg$gamma_slopes
  
  if (is.null(SigmaX)) {
    SigmaX <- make_SigmaX(ar = 0.3)
  }
  
  p <- length(beta)
  if (ncol(SigmaX) != p) {
    stop("SigmaX dimension does not match length(beta).")
  }
  
  # --- generate X ~ N(0, SigmaX) ---
  Lx <- chol(SigmaX)
  Zx <- matrix(rnorm(n * p), n, p)
  X  <- Zx %*% Lx   # n x p
  
  # --- build Z = (1, X) for selection ---
  Z <- cbind(1, X)
  colnames(X) <- paste0("x", seq_len(p))
  colnames(Z) <- c("intercept", paste0("x", seq_len(p)))
  
  # --- build gamma = (gamma0, gamma_slopes) ---
  if (is.null(scen$gamma0) || is.na(scen$gamma0)) {
    stop("scen$gamma0 is missing - did you run add_gamma0_to_grid()?")
  }
  gamma0 <- scen$gamma0
  gamma  <- c(gamma0, gamma_slopes)
  
  # --- draw (u, v) according to error family ---
  efam <- as.character(scen$err_family)
  
  errs <- switch(
    efam,
    "normal" = draw_uv_normal(n, rho),
    "t3"     = draw_uv_t(n, rho, df = 3),
    "t5"     = draw_uv_t(n, rho, df = 5),
    "mixture" = draw_uv_mixture(
      n = n,
      rho = rho,
      base_family = "normal",           # v ~ normal, u gets mixture shift
      df = NULL,
      mix_pi = scen$mix_pi,
      mix_shift_u = scen$mix_shift_u
    ),
    stop("Unknown err_family in simulate_one_dataset: ", efam)
  )
  
  u <- errs$u
  v <- errs$v
  
  # --- latent outcome and selection index ---
  y_star <- as.numeric(X %*% beta + u)
  s_star <- as.numeric(Z %*% gamma + v)
  
  # --- selection and observed Y ---
  S <- as.integer(s_star > 0)
  Y <- ifelse(S == 1L, y_star, NA_real_)
  
  # --- return dataset ---
  list(
    Y = Y,
    S = S,
    X = X,
    Z = Z,
    y_star = y_star,
    s_star = s_star
  )
}
