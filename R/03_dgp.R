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

# --- helper: draw bivariate normal errors (u, v) with corr = rho ---
draw_uv_normal <- function(n, rho) {
  Sigma <- matrix(c(1, rho,
                    rho, 1), 2, 2)
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * 2), n, 2)
  E <- Z %*% L
  list(u = E[, 1], v = E[, 2])
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

# --- helper: draw "logistic" errors using a Gaussian copula ---
# We want marginals ~ standard logistic (mean 0, var 1) and correlation ~ rho.
# Construction:
#   1) draw bivariate normal with corr = rho
#   2) map to U(0,1) via Phi
#   3) map to logistic via qlogis
#   4) rescale to variance 1 (standard logistic has var = pi^2 / 3)
draw_uv_logistic <- function(n, rho) {
  # step 1: correlated normals
  Sigma <- matrix(c(1, rho,
                    rho, 1), 2, 2)
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * 2), n, 2)
  N2 <- Z %*% L
  
  # step 2: map to uniforms
  U  <- pnorm(N2)
  
  # step 3: map to logistic(0, 1) (var = pi^2 / 3)
  Lraw <- qlogis(U)
  
  # step 4: rescale to variance ~ 1
  scale <- sqrt(3) / pi
  Lstd  <- Lraw * scale
  
  list(u = Lstd[, 1], v = Lstd[, 2])
}

simulate_one_dataset <- function(scen, cfg, SigmaX = NULL) {
  # scen: one-row data.frame with n, rho, err_family, df, gamma0, (sel_model not used here)
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
    "normal"   = draw_uv_normal(n, rho),
    "t3"       = draw_uv_t(n, rho, df = 3),
    "t5"       = draw_uv_t(n, rho, df = 5),
    "logistic" = draw_uv_logistic(n, rho),
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
    Y      = Y,
    S      = S,
    X      = X,
    Z      = Z,
    y_star = y_star,
    s_star = s_star
  )
}
