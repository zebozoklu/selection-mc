# R/05_mc_engine.R
# Monte Carlo engine
# Purpose: repeat simulation + estimation R times, for each scenario
# Inputs: scenario table, number of reps R, estimator functions
# Outputs:
#   - raw_beta: estimates per rep (R x K) per estimator
#   - raw_se:   robust SEs per rep (R x K) per estimator
#   - aggregated summaries per scenario/estimator
# Sanity checks: reproducible with seed; progress logging; stable coef naming

run_mc_one_scenario <- function(scen, cfg, estimators,
                                R = cfg$R_debug,
                                seed = cfg$seed_mc) {
  
  set.seed(seed)
  
  p <- length(cfg$beta)       # number of X slopes
  K <- p + 1                  # intercept + slopes
  
  # Canonical coefficient names (match your estimators: "intercept", "x1", "x2", ...)
  coef_names <- c("intercept", paste0("x", seq_len(p)))
  
  # Coerce estimator output into canonical order/names
  coerce_to_canonical <- function(v, what = "beta") {
    if (is.null(v)) stop("Estimator returned NULL ", what, ".")
    if (length(v) != K) {
      stop("Estimator returned ", what, " length ", length(v),
           " but expected ", K, ".")
    }
    
    # If unnamed, assume it's already in canonical order
    if (is.null(names(v)) || any(names(v) == "")) {
      names(v) <- coef_names
      return(v)
    }
    
    # Normalize common intercept name variants -> "intercept"
    nm <- names(v)
    nm[nm %in% c("(Intercept)", "Intercept", "CONST", "const", "Constant", "constant")] <- "intercept"
    names(v) <- nm
    
    # Require all canonical names and reorder
    if (!all(coef_names %in% names(v))) {
      stop("Estimator returned ", what, " with names {",
           paste(names(v), collapse = ", "),
           "} but expected at least {",
           paste(coef_names, collapse = ", "),
           "}.")
    }
    
    v[coef_names]
  }
  
  # Prepare containers with fixed colnames (critical for density plots)
  raw_beta <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K, dimnames = list(NULL, coef_names))
  })
  raw_se <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K, dimnames = list(NULL, coef_names))
  })
  
  # Optional diagnostic
  sel_rate <- rep(NA_real_, R)
  
  # Main loop
  for (r in seq_len(R)) {
    dat <- simulate_one_dataset(scen, cfg)
    
    sel_rate[r] <- mean(dat$S)
    
    for (nm in names(estimators)) {
      out <- estimators[[nm]](dat)
      
      if (!is.list(out) || is.null(out$beta) || is.null(out$se)) {
        stop("Estimator ", nm, " must return list(beta=..., se=...).")
      }
      
      b_hat  <- coerce_to_canonical(out$beta, what = "beta")
      se_hat <- coerce_to_canonical(out$se,   what = "se")
      
      raw_beta[[nm]][r, ] <- b_hat
      raw_se[[nm]][r, ]   <- se_hat
    }
  }
  
  # True coefficients (assume true intercept = 0), named to match canonical scheme
  beta_true <- setNames(c(0, cfg$beta), coef_names)
  
  # Summarize
  summary <- lapply(names(raw_beta), function(nm) {
    summarize_mc_estimator(
      beta_hat_mat = raw_beta[[nm]],
      se_hat_mat   = raw_se[[nm]],
      beta_true    = beta_true,
      # cfg$target_beta_index is index among slopes: x1=1, x2=2, ...
      target_index = 1 + cfg$target_beta_index
    )
  })
  names(summary) <- names(raw_beta)
  
  list(
    raw_beta = raw_beta,
    raw_se   = raw_se,
    sel_rate = sel_rate,
    summary  = summary,
    scen     = scen
  )
}


