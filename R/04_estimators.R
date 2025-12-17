# Estimator layer
# Purpose: implement each estimator mapping (Y,S,X,Z) -> beta_hat
# Estimators:
#   - selected-sample OLS
#   - zero-imputation OLS
#   - Heckman two-step variants (probit/logit/LPM first stage)
# Inputs: simulated dataset list
# Outputs: standardized beta_hat vector (intercept + slopes)
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

# --- helper: common second stage given a selection index z_index_hat ---
# We always build a "Heckman-style" control function using the normal IMR:
#   lambda_i = phi(z_i) / Phi(z_i)
# even if the first stage is logit or LPM -> deliberately mis-specified.
heckman_second_stage <- function(dat, z_index_hat) {
  # Normal-based IMR from supplied index
  p_hat <- pnorm(z_index_hat)
  # avoid division by zero / Inf in extreme cases
  eps <- 1e-8
  p_hat <- pmin(pmax(p_hat, eps), 1 - eps)
  lambda <- dnorm(z_index_hat) / p_hat
  
  # restrict to selected sample
  idx <- which(dat$S == 1L)
  if (length(idx) == 0L) stop("No selected observations for Heckman two-step.")
  
  Ys    <- dat$Y[idx]
  Xs    <- dat$X[idx, , drop = FALSE]
  lam_s <- lambda[idx]
  
  X_design <- cbind(1, Xs, lam_s)
  
  b_full <- ols_fit(Ys, X_design)
  
  # first 1 + ncol(X) entries correspond to intercept + beta's
  p <- ncol(dat$X)
  b <- b_full[1:(1 + p)]
  names(b) <- c("intercept", colnames(dat$X))
  
  b
}

# --- 3a) Heckman two-step with PROBIT first stage (the classical one) ---
est_heckman_probit <- function(dat) {
  # 1) Probit selection: S ~ Z
  sel_df <- data.frame(S = dat$S, dat$Z)
  # dat$Z already includes intercept column; fit without an extra intercept
  sel_fit <- glm(S ~ . - 1, data = sel_df,
                 family = binomial(link = "probit"))
  
  # linear index (on probit scale)
  z_index_hat <- as.numeric(dat$Z %*% coef(sel_fit))
  
  # 2) Common second stage
  heckman_second_stage(dat, z_index_hat)
}

# keep old name for backward compatibility
est_heckman_2step <- est_heckman_probit

# --- 3b) "Heckman-logit": logit first stage + normal IMR control function ---
est_heckman_logit <- function(dat) {
  sel_df <- data.frame(S = dat$S, dat$Z)
  
  sel_fit <- glm(S ~ . - 1, data = sel_df,
                 family = binomial(link = "logit"))
  
  # linear predictor on logit scale
  z_index_hat <- as.numeric(dat$Z %*% coef(sel_fit))
  # Alternatively: z_index_hat <- predict(sel_fit, type = "link")
  
  heckman_second_stage(dat, z_index_hat)
}

# --- 3c) "Heckman-LPM": linear probability first stage + normal IMR ---
est_heckman_lpm <- function(dat) {
  sel_df <- data.frame(S = dat$S, dat$Z)
  
  sel_fit <- lm(S ~ . - 1, data = sel_df)
  
  z_index_hat <- as.numeric(dat$Z %*% coef(sel_fit))
  
  heckman_second_stage(dat, z_index_hat)
}



