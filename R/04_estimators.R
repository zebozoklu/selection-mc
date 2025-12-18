# Estimator layer
# Purpose: implement each estimator mapping (Y,S,X,Z) -> beta_hat (+ se)
# Estimators:
#   - selected-sample OLS
#   - zero-imputation OLS
#   - Heckman two-step variants (probit/logit/LPM first stage)
# Inputs: simulated dataset list + cfg
# Outputs: list(beta = ..., se = ...)  (intercept + slopes)
# Sanity checks: correct length; no silent NA

# --- helper: OLS with either homoskedastic or HC1 robust SEs using lm.fit ---
ols_fit_with_se <- function(y, X, se_type = c("HC1", "homosked")) {
  se_type <- match.arg(se_type)
  
  fit <- lm.fit(x = X, y = y)
  b   <- coef(fit)
  
  n <- NROW(X)
  k <- NCOL(X)
  
  # residuals
  e <- as.numeric(y - X %*% b)
  
  XtX_inv <- solve(crossprod(X))
  
  if (se_type == "homosked") {
    sigma2 <- sum(e^2) / (n - k)
    vcov_b <- sigma2 * XtX_inv
  } else {
    # HC1 sandwich: (n/(n-k)) * (X'X)^(-1) X' diag(e^2) X (X'X)^(-1)
    meat   <- crossprod(X, X * (e^2))
    vcov_b <- (n / (n - k)) * XtX_inv %*% meat %*% XtX_inv
  }
  
  se <- sqrt(diag(vcov_b))
  
  list(beta = as.numeric(b), se = as.numeric(se))
}

# --- 1) Selected-sample OLS (ignore selection) ---
est_selected_ols <- function(dat, cfg) {
  idx <- which(dat$S == 1L)
  if (length(idx) == 0L) stop("No selected observations for selected-OLS.")
  
  Xs <- dat$X[idx, , drop = FALSE]
  Ys <- dat$Y[idx]
  
  X_design <- cbind(1, Xs)  # add intercept
  
  fit <- ols_fit_with_se(Ys, X_design, se_type = cfg$se_type)
  b   <- fit$beta
  se  <- fit$se
  
  names(b)  <- c("intercept", colnames(dat$X))
  names(se) <- c("intercept", colnames(dat$X))
  
  list(beta = b, se = se)
}

# --- 2) Zero-imputation OLS ---
est_zero_impute_ols <- function(dat, cfg) {
  Y0 <- dat$Y
  Y0[is.na(Y0)] <- 0
  
  X_design <- cbind(1, dat$X)
  
  fit <- ols_fit_with_se(Y0, X_design, se_type = cfg$se_type)
  b   <- fit$beta
  se  <- fit$se
  
  names(b)  <- c("intercept", colnames(dat$X))
  names(se) <- c("intercept", colnames(dat$X))
  
  list(beta = b, se = se)
}

# --- helper: second stage using normal IMR, built from predicted probabilities ---
# We compute:
#   z_tilde = Phi^{-1}(p_hat)
#   lambda  = phi(z_tilde) / p_hat
# This avoids mixing scales (logit/LPM indexes are not probit indexes).
heckman_second_stage_from_phat <- function(dat, p_hat, cfg) {
  eps <- cfg$eps_prob %||% 1e-8
  p_hat <- pmin(pmax(p_hat, eps), 1 - eps)
  
  z_tilde <- qnorm(p_hat)
  lambda  <- dnorm(z_tilde) / p_hat
  
  idx <- which(dat$S == 1L)
  if (length(idx) == 0L) stop("No selected observations for Heckman two-step.")
  
  Ys    <- dat$Y[idx]
  Xs    <- dat$X[idx, , drop = FALSE]
  lam_s <- lambda[idx]
  
  X_design <- cbind(1, Xs, lam_s)
  
  fit <- ols_fit_with_se(Ys, X_design, se_type = cfg$se_type)
  b_full  <- fit$beta
  se_full <- fit$se
  
  p <- ncol(dat$X)
  keep <- 1:(1 + p)
  
  b  <- b_full[keep]
  se <- se_full[keep]
  
  names(b)  <- c("intercept", colnames(dat$X))
  names(se) <- c("intercept", colnames(dat$X))
  
  list(beta = b, se = se)
}

# --- 3a) Heckman two-step with PROBIT first stage ---
est_heckman_probit <- function(dat, cfg) {
  sel_df <- data.frame(S = dat$S, dat$Z)
  
  sel_fit <- glm(
    S ~ . - 1,
    data   = sel_df,
    family = binomial(link = "probit")
  )
  
  p_hat <- fitted(sel_fit)  # predicted selection prob
  heckman_second_stage_from_phat(dat, p_hat, cfg)
}

# keep old name for backward compatibility
est_heckman_2step <- est_heckman_probit

# --- 3b) Heckman "logit": logit first stage + normal IMR from p_hat ---
est_heckman_logit <- function(dat, cfg) {
  sel_df <- data.frame(S = dat$S, dat$Z)
  
  sel_fit <- glm(
    S ~ . - 1,
    data   = sel_df,
    family = binomial(link = "logit")
  )
  
  p_hat <- fitted(sel_fit)
  heckman_second_stage_from_phat(dat, p_hat, cfg)
}

# --- 3c) Heckman "LPM": linear probability first stage + normal IMR from p_hat ---
est_heckman_lpm <- function(dat, cfg) {
  sel_df <- data.frame(S = dat$S, dat$Z)
  
  sel_fit <- lm(S ~ . - 1, data = sel_df)
  
  p_hat <- as.numeric(dat$Z %*% coef(sel_fit))
  eps <- cfg$eps_prob %||% 1e-8
  p_hat <- pmin(pmax(p_hat, eps), 1 - eps)
  
  heckman_second_stage_from_phat(dat, p_hat, cfg)
}



