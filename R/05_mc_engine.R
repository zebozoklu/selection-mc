# R/05_mc_engine.R

run_mc_one_scenario <- function(scen, cfg, estimators,
                                R = cfg$R_debug,
                                seed = cfg$seed_mc) {
  
  set.seed(seed)
  
  p <- length(cfg$beta)       # number of X slopes
  K <- p + 1                  # intercept + slopes
  
  # Canonical coefficient names used everywhere
  coef_names <- c("(Intercept)", paste0("x", seq_len(p)))
  
  # prepare containers WITH colnames so density plots always work
  raw_beta <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K, dimnames = list(NULL, coef_names))
  })
  raw_se <- lapply(estimators, function(f) {
    matrix(NA_real_, nrow = R, ncol = K, dimnames = list(NULL, coef_names))
  })
  
  # optional diagnostic
  sel_rate <- rep(NA_real_, R)
  
  # helper: coerce estimator output into canonical order/names
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
    
    # If named: require all canonical names and reorder
    if (!all(coef_names %in% names(v))) {
      stop("Estimator returned ", what, " with names {",
           paste(names(v), collapse = ", "),
           "} but expected at least {",
           paste(coef_names, collapse = ", "),
           "}.")
    }
    
    v[coef_names]
  }
  
  # main loop
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
  
  # true coefficients (assume true intercept = 0) — NAME THEM
  beta_true <- setNames(c(0, cfg$beta), coef_names)
  
  # summarize
  summary <- lapply(names(raw_beta), function(nm) {
    summarize_mc_estimator(
      beta_hat_mat = raw_beta[[nm]],
      se_hat_mat   = raw_se[[nm]],
      beta_true    = beta_true,
      # if cfg$target_beta_index refers to x1=1, x2=2, ... then +1 accounts for intercept
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


