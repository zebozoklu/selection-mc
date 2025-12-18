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

# --- helper: bare-bones OLS using lm.fit ---
ols_fit <- function(y, X) {
  fit <- lm.fit(x = X, y = y)
  out <- coef(fit)
  # ensure it's a plain numeric, not named weirdly
  as.numeric(out)
}

# --- 1) Selected-sample OLS (ignore selection) ---
est_selected_ols <- function(dat) {
  idx <- which(dat$S == 1L)
  if (length(idx) == 0L) stop("No selected observations for selected-OLS.")
  
  Xs <- dat$X[idx, , drop = FALSE]
  Ys <- dat$Y[idx]
  
  X_design <- cbind(1, Xs)       # add intercept
  b <- ols_fit(Ys, X_design)
  
  names(b) <- c("intercept", colnames(dat$X))
  b
}

# --- 2) Zero-imputation OLS ---
est_zero_impute_ols <- function(dat) {
  Y0 <- dat$Y
  # replace NA outcomes with zero
  Y0[is.na(Y0)] <- 0
  
  X_design <- cbind(1, dat$X)
  b <- ols_fit(Y0, X_design)
  
  names(b) <- c("intercept", colnames(dat$X))
  b
}

# --- 3) Heckman two-step ---
est_heckman_2step <- function(dat) {
  # 1) Probit selection: S ~ Z
  sel_df <- data.frame(S = dat$S, dat$Z)
  # dat$Z already includes intercept column; fit without an extra intercept
  sel_fit <- glm(S ~ . - 1, data = sel_df,
                 family = binomial(link = "probit"))
  
  z_index_hat <- as.numeric(dat$Z %*% coef(sel_fit))
  
  # Inverse Mills ratio
  p_hat <- pnorm(z_index_hat)
  lambda <- dnorm(z_index_hat) / p_hat
  
  # restrict to selected sample
  idx <- which(dat$S == 1L)
  if (length(idx) == 0L) stop("No selected observations for Heckman two-step.")
  
  Ys  <- dat$Y[idx]
  Xs  <- dat$X[idx, , drop = FALSE]
  lam_s <- lambda[idx]
  
  X_design <- cbind(1, Xs, lam_s)
  
  b_full <- ols_fit(Ys, X_design)
  
  # first 1 + ncol(X) entries correspond to intercept + beta's
  p <- ncol(dat$X)
  b <- b_full[1:(1 + p)]
  names(b) <- c("intercept", colnames(dat$X))
  
  b
}
