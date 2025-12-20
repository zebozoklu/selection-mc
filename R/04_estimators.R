# Estimator layer
# Purpose: implement each estimator mapping (Y,S,X,Z) -> beta_hat + robust SE
# Estimators:
#   - selected-sample OLS
#   - zero-imputation OLS
#   - Heckman two-step (probit + IMR) with stacked-sandwich SE
# Outputs: list(beta = ..., se = ..., aux = optional)

# ---------- small utilities ----------

safe_ginv <- function(A) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Need MASS for ginv fallback when a matrix is singular.")
  }
  MASS::ginv(A)
}

safe_solve <- function(A) {
  tryCatch(solve(A), error = function(e) safe_ginv(A))
}

# Stable inverse Mills ratio for probit:
# lambda(z) = phi(z) / Phi(z) computed as exp(log phi - log Phi)
imr_probit <- function(z, eps = 1e-12) {
  # clamp Phi(z) in (eps, 1-eps) to avoid log(0) and division blowups
  p <- pnorm(z)
  p <- pmin(pmax(p, eps), 1 - eps)
  
  log_phi <- dnorm(z, log = TRUE)
  log_Phi <- pnorm(z, log.p = TRUE)  # log Phi(z)
  
  # exp(log_phi - log_Phi) is stable for very negative z
  lam <- exp(log_phi - log_Phi)
  
  # additional cap just in case (prevents catastrophic inf from eps)
  lam <- pmin(lam, 1 / eps)
  lam
}

# ---------- helper: OLS + robust (HC1) SE using lm.fit ----------

ols_fit_robust <- function(y, X, hc = c("HC1", "HC0")) {
  hc <- match.arg(hc)
  
  fit <- lm.fit(x = X, y = y)
  b <- as.numeric(coef(fit))
  
  if (anyNA(b)) {
    stop("OLS produced NA coefficients (rank deficiency / collinearity).")
  }
  
  e <- as.numeric(fit$residuals)
  n <- NROW(X)
  k <- NCOL(X)
  
  XtX <- crossprod(X)
  
  # Invert X'X (fallback to generalized inverse if needed)
  XtX_inv <- tryCatch(
    solve(XtX),
    error = function(err) safe_ginv(XtX)
  )
  
  # meat = X' diag(e^2) X
  meat <- crossprod(X, X * (e^2))
  
  V <- XtX_inv %*% meat %*% XtX_inv
  
  if (hc == "HC1") {
    V <- (n / (n - k)) * V
  }
  
  se <- sqrt(pmax(diag(V), 0))
  
  list(beta = b, se = as.numeric(se), vcov = V)
}

# ---------- 1) Selected-sample OLS (ignore selection) ----------

est_selected_ols <- function(dat) {
  idx <- which(dat$S == 1L)
  if (length(idx) == 0L) stop("No selected observations for selected-OLS.")
  
  Xs <- dat$X[idx, , drop = FALSE]
  Ys <- dat$Y[idx]
  
  X_design <- cbind(1, Xs)  # add intercept
  
  fit <- ols_fit_robust(Ys, X_design, hc = "HC1")
  b <- fit$beta
  se <- fit$se
  
  nm <- c("intercept", colnames(dat$X))
  names(b) <- nm
  names(se) <- nm
  
  list(beta = b, se = se)
}

# ---------- 2) Zero-imputation OLS ----------

est_zero_impute_ols <- function(dat) {
  Y0 <- dat$Y
  Y0[is.na(Y0)] <- 0
  
  X_design <- cbind(1, dat$X)
  
  fit <- ols_fit_robust(Y0, X_design, hc = "HC1")
  b <- fit$beta
  se <- fit$se
  
  nm <- c("intercept", colnames(dat$X))
  names(b) <- nm
  names(se) <- nm
  
  list(beta = b, se = se)
}

# ---------- 3) Heckman two-step (probit + IMR) with stacked sandwich SE ----------
# This implements a stacked-moment (Murphy–Topel / GMM-style) sandwich:
#   m1_i(gamma) = probit score moment
#   m2_i(gamma,theta2) = I(S=1) * w_i(gamma) * u_i(gamma,theta2)
# and V = (G^{-1} Omega G^{-1'}) / n using numerical derivatives for G11, G21.

est_heckman_2step <- function(dat, eps = 1e-10, fd_step = 1e-6) {
  Z <- as.matrix(dat$Z)            # includes intercept already
  S <- as.numeric(dat$S)
  n <- NROW(Z)
  
  idx <- which(S == 1)
  if (length(idx) == 0L) stop("No selected observations for Heckman two-step.")
  
  Xs <- dat$X[idx, , drop = FALSE]
  Ys <- dat$Y[idx]
  
  # ----- Step 1: probit selection -----
  # Use glm.fit directly to avoid data.frame column-name quirks.
  sel_fit <- glm.fit(
    x = Z,
    y = S,
    family = binomial(link = "probit")
  )
  gamma_hat <- as.numeric(coef(sel_fit))
  if (anyNA(gamma_hat)) stop("Probit produced NA coefficients (separation / collinearity).")
  
  eta <- as.numeric(Z %*% gamma_hat)
  
  # Clamp p only for score factors; lambda computed stably via logs
  p <- pnorm(eta)
  p <- pmin(pmax(p, eps), 1 - eps)
  phi <- dnorm(eta)
  
  lambda <- imr_probit(eta, eps = eps)
  
  # ----- Step 2: OLS on selected with IMR -----
  W <- cbind(1, Xs, lambda[idx])          # [1, X, lambda]
  k2 <- NCOL(W)
  
  fit2 <- lm.fit(x = W, y = Ys)
  theta2_hat <- as.numeric(coef(fit2))   # includes lambda coefficient last
  if (anyNA(theta2_hat)) stop("Second-stage OLS produced NA coefficients (collinearity).")
  
  u_sel <- as.numeric(fit2$residuals)
  
  # Build w_all and u_all as n-length objects (0 for nonselected) for moment stacking
  w_all <- matrix(0, n, k2)
  w_all[idx, ] <- W
  u_all <- rep(0, n)
  u_all[idx] <- u_sel
  
  # ----- Moment vectors -----
  # m1_i (probit score moment): Z_i * [ (S - p) * phi / (p(1-p)) ]
  score_fac <- (S - p) * phi / (p * (1 - p))
  m1 <- Z * score_fac
  k1 <- NCOL(m1)
  
  # m2_i (OLS moment on selected): I(S=1) * w_i * u_i
  m2 <- w_all * u_all
  
  m <- cbind(m1, m2)
  Omega_hat <- crossprod(m) / n
  
  # ----- Jacobian blocks G11, G21 numerically; G22 analytically -----
  
  mean_score <- function(gamma) {
    eta_g <- as.numeric(Z %*% gamma)
    p_g <- pnorm(eta_g)
    p_g <- pmin(pmax(p_g, eps), 1 - eps)
    phi_g <- dnorm(eta_g)
    fac_g <- (S - p_g) * phi_g / (p_g * (1 - p_g))
    colMeans(Z * fac_g)
  }
  
  mean_m2 <- function(gamma) {
    eta_g <- as.numeric(Z %*% gamma)
    # lambda needs stable computation
    lam_g <- imr_probit(eta_g, eps = eps)
    
    W_g <- cbind(1, Xs, lam_g[idx])
    u_g <- as.numeric(Ys - W_g %*% theta2_hat)
    
    w_all_g <- matrix(0, n, k2)
    w_all_g[idx, ] <- W_g
    u_all_g <- rep(0, n)
    u_all_g[idx] <- u_g
    
    colMeans(w_all_g * u_all_g)
  }
  
  # choose scale-aware finite-difference step
  step_j <- function(xj) fd_step * max(1, abs(xj))
  
  G11 <- matrix(0, k1, k1)
  G21 <- matrix(0, k2, k1)
  
  for (j in seq_len(k1)) {
    h <- step_j(gamma_hat[j])
    e_j <- rep(0, k1); e_j[j] <- 1
    
    gp <- gamma_hat + h * e_j
    gm <- gamma_hat - h * e_j
    
    G11[, j] <- (mean_score(gp) - mean_score(gm)) / (2 * h)
    G21[, j] <- (mean_m2(gp)    - mean_m2(gm))    / (2 * h)
  }
  
  # G22: derivative wrt theta2 of E[I(S=1) w (Y - w'theta2)] = -E[I(S=1) w w']
  G22 <- - crossprod(W) / n
  
  # Assemble G
  G_hat <- matrix(0, k1 + k2, k1 + k2)
  G_hat[1:k1, 1:k1] <- G11
  G_hat[(k1 + 1):(k1 + k2), 1:k1] <- G21
  G_hat[(k1 + 1):(k1 + k2), (k1 + 1):(k1 + k2)] <- G22
  
  # V = (G^{-1} Omega G^{-1'}) / n
  Ginv <- safe_solve(G_hat)
  V_hat <- (Ginv %*% Omega_hat %*% t(Ginv)) / n
  
  # Outcome-parameter block (theta2)
  V_theta2 <- V_hat[(k1 + 1):(k1 + k2), (k1 + 1):(k1 + k2)]
  se_theta2 <- sqrt(pmax(diag(V_theta2), 0))
  
  # Return only intercept + slopes (drop lambda coefficient), consistent with your API
  pX <- ncol(dat$X)
  b_full  <- theta2_hat
  se_full <- se_theta2
  
  b  <- b_full[1:(1 + pX)]
  se <- se_full[1:(1 + pX)]
  
  nm <- c("intercept", colnames(dat$X))
  names(b)  <- nm
  names(se) <- nm
  
  list(
    beta = b,
    se   = se,
    aux  = list(
      se_type     = "stacked_sandwich",
      gamma_hat   = gamma_hat,
      lambda_coef = b_full[1 + pX + 1],
      lambda_se   = se_full[1 + pX + 1]
    )
  )
}


