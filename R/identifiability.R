# ============================================================
# R/identifiability.R
#
# Orchestration for the A1 identifiability analysis:
#   - single_recovery()    : one simulate-then-estimate experiment
#   - recovery_study()     : N-replication study with summary statistics
#   - plot_recovery()      : dot-plot of bias and RMSE per parameter
#   - plot_profiles()      : 2x3 grid of profile likelihood curves
#
# Requires: model_functions.R and mle.R sourced first.
# ============================================================

# ---- Seeds for reproducibility ----
# The master seed produces all per-replication seeds deterministically.
RECOVERY_MASTER_SEED <- 2026L

# data.table is used for fast row-binding in summary construction.
if (!requireNamespace("data.table", quietly = TRUE))
  stop("Package 'data.table' required. Install with install.packages('data.table').")

# ---- Canonical recovery protocol ----
# Used uniformly across ALL Monte Carlo studies (conditional, marginal,
# duration, sensitivity) to ensure a single principled, physiologically
# realistic design.  The double-logistic transient stressor:
#   - is smooth and differentiable everywhere
#   - creates sufficient non-stationarity to decouple (σ_p, σ_s) from (μ₀, ρ)
#   - has known c_p, c_s treated as fixed constants during inference
#   - makes θ_true = θ* so recovery metrics measure estimator precision only

CANONICAL_INPUT_FN <- make_double_logistic(t_on = 300, t_off = 420, k = 0.1)
CANONICAL_DURATION <- 720L   # 300 s baseline + 120 s stress + 300 s recovery

# ---- Single recovery experiment ----
#
# Simulates one dataset from true_params, runs the conditional MLE,
# and returns true values, estimates, and SEs side-by-side.
#
# NOTE: this is *conditional* identifiability, meaning parameters are estimated
# given the true state trajectory. Full marginal recovery (from spike
# times alone) requires UKF

single_recovery <- function(true_params,
                            duration = CANONICAL_DURATION, dt = 0.005,
                            input_fn = CANONICAL_INPUT_FN,
                            seed = NULL) {
  sim_res <- sim_sde_ig(duration, dt, true_params, input_fn, seed = seed)
  # Pass input_fn so the conditional MLE removes the known input contribution
  # from OU residuals.  c_p and c_s are taken from true_params (fixed, not estimated).
  mle     <- full_conditional_mle(sim_res, input_fn = input_fn)

  fp           <- true_params$free
  use_coupled  <- (abs(fp$a_ps) + abs(fp$a_sp)) > 1e-12
  true_vec <- c(a_p = fp$a_p, a_s = fp$a_s,
                sigma_p = fp$sigma_p, sigma_s = fp$sigma_s,
                mu0 = fp$mu_0, rho = fp$rho)
  hat_vec  <- c(a_p = mle$a_p,  a_s = mle$a_s,
                sigma_p = mle$sigma_p, sigma_s = mle$sigma_s,
                mu0 = mle$mu0, rho = mle$rho)
  se_vec   <- c(a_p = mle$a_p_se, a_s = mle$a_s_se,
                sigma_p = mle$sigma_p_se, sigma_s = mle$sigma_s_se,
                mu0 = mle$mu0_se, rho = mle$rho_se)
  if (use_coupled) {
    true_vec <- c(true_vec, a_ps = fp$a_ps, a_sp = fp$a_sp)
    hat_vec  <- c(hat_vec,  a_ps = mle$a_ps, a_sp = mle$a_sp)
    se_vec   <- c(se_vec,   a_ps = mle$a_ps_se, a_sp = mle$a_sp_se)
  }

  list(true = true_vec, hat = hat_vec, se = se_vec,
       n_beats = mle$n_beats, ll = mle$ll_total)
}

# ---- N-replication recovery study ----
#
# Parameters
#   N            : number of Monte Carlo replications
#   true_params  : make_model_params() output
#   duration     : simulation length in seconds
#   dt           : time step
#   input_fn     : exogenous input; must be non-trivial to identify c_p, c_s
#   use_parallel : use parallel::mclapply on Unix; falls back to lapply

recovery_study <- function(N = 200L,
                           true_params,
                           duration = CANONICAL_DURATION,
                           dt = 0.005,
                           input_fn = CANONICAL_INPUT_FN,
                           use_parallel = TRUE) {
  set.seed(RECOVERY_MASTER_SEED)
  seeds <- sample.int(1e6L, N)

  run_one <- function(i) {
    tryCatch(
      single_recovery(true_params, duration, dt, input_fn, seed = seeds[i]),
      error = function(e) NULL
    )
  }

  if (use_parallel && .Platform$OS.type == "unix" &&
      requireNamespace("parallel", quietly = TRUE)) {
    n_cores  <- max(1L, parallel::detectCores() - 1L)
    results  <- parallel::mclapply(seq_len(N), run_one, mc.cores = n_cores,
                                   mc.set.seed = FALSE)
  } else {
    results  <- lapply(seq_len(N), run_one)
  }
  results <- Filter(Negate(is.null), results)

  use_coupled <- (abs(true_params$free$a_ps) + abs(true_params$free$a_sp)) > 1e-12
  params <- c("a_p", "a_s", "sigma_p", "sigma_s", "mu0", "rho")
  labels <- c(expression(a[p]), expression(a[s]),
              expression(sigma[p]), expression(sigma[s]),
              expression(mu[0]), expression(rho))
  if (use_coupled) {
    params <- c(params, "a_ps", "a_sp")
    labels <- c(labels, expression(a[ps]), expression(a[sp]))
  }

  summary_df <- data.table::rbindlist(lapply(seq_along(params), function(j) {
    p      <- params[j]
    true_v <- sapply(results, function(r) r$true[p])
    hat_v  <- sapply(results, function(r) r$hat[p])
    se_v   <- sapply(results, function(r) r$se[p])

    # Separate validity: bias/RMSE only require finite estimates;
    # coverage additionally requires finite, positive SEs.
    ok_val <- is.finite(hat_v) & is.finite(true_v)
    ok_cov <- ok_val & is.finite(se_v) & se_v > 0

    n_ok  <- sum(ok_val)
    n_cov <- sum(ok_cov)
    if (n_ok == 0L) {
      return(data.frame(
        parameter = p, true_value = NA_real_, mean_hat = NA_real_,
        bias = NA_real_, rel_bias_pct = NA_real_,
        rmse = NA_real_, rmse_rel_pct = NA_real_,
        coverage_95 = NA_real_, n_valid = 0L,
        stringsAsFactors = FALSE))
    }
    t_bar <- mean(true_v[ok_val])
    h_bar <- mean(hat_v[ok_val])
    bias  <- h_bar - t_bar
    rmse  <- sqrt(mean((hat_v[ok_val] - true_v[ok_val])^2))

    # 95% Wald CI coverage — only computable when SEs are finite and positive
    cover <- if (n_cov > 0L) {
      ci_lo <- hat_v[ok_cov] - 1.96 * se_v[ok_cov]
      ci_hi <- hat_v[ok_cov] + 1.96 * se_v[ok_cov]
      mean(ci_lo <= true_v[ok_cov] & true_v[ok_cov] <= ci_hi)
    } else NA_real_

    data.frame(
      parameter     = p,
      true_value    = t_bar,
      mean_hat      = h_bar,
      bias          = bias,
      rel_bias_pct  = 100 * bias / abs(t_bar),
      rmse          = rmse,
      rmse_rel_pct  = 100 * rmse / abs(t_bar),
      coverage_95   = round(100 * cover, 1),
      n_valid       = n_ok,       # replications contributing to bias/RMSE
      n_cov_valid   = n_cov,      # replications contributing to coverage (SE finite & > 0)
      stringsAsFactors = FALSE
    )
  }))

  summary_df <- as.data.frame(summary_df)

  list(
    results    = results,
    summary    = summary_df,
    N_total    = N,
    N_valid    = length(results),
    true_params = true_params
  )
}

# ---- Recovery dot-plot ----
#
# Displays relative bias (center) ± relative RMSE (whiskers) per parameter.
# Inputs: the summary data.frame from recovery_study().

plot_recovery <- function(rec_summary,
                          main = "Conditional MLE Recovery Study") {
  df  <- rec_summary
  par_labels <- c(
    a_p     = expression(a[p]~"(Hz)"),
    a_s     = expression(a[s]~"(Hz)"),
    sigma_p = expression(sigma[p]),
    sigma_s = expression(sigma[s]),
    mu0     = expression(mu[0]),
    rho     = expression(rho~"(CV)"),
    a_ps    = expression(a[ps]),
    a_sp    = expression(a[sp])
  )

  np   <- nrow(df)
  base_cols <- c(
    a_p     = "#0072B2",
    a_s     = "#56B4E9",
    sigma_p = "#D55E00",
    sigma_s = "#E69F00",
    mu0     = "#2C5F2E",
    rho     = "#2C5F2E",
    a_ps    = "#999999",
    a_sp    = "#999999"
  )
  cols <- unname(base_cols[df$parameter])
  cols[is.na(cols)] <- "black"

  par(mar = c(4.5, 6, 3, 5))
  xlim_vals <- c(df$rel_bias_pct - df$rmse_rel_pct,
                 df$rel_bias_pct + df$rmse_rel_pct)
  xlim_vals <- xlim_vals[is.finite(xlim_vals)]
  if (length(xlim_vals) == 0L) xlim_vals <- c(-10, 10)
  plot(df$rel_bias_pct, seq_len(np), pch = 19, cex = 1.5,
       col  = cols,
       xlim = range(xlim_vals) * 1.15,
       ylim = c(0.5, np + 0.5),
       xlab = "Relative error (%)",
       ylab = "",
       yaxt = "n",
       main = main,
       frame.plot = FALSE)
  abline(v = 0, lty = 2, col = "gray50")
  for (i in seq_len(np)) {
    segments(df$rel_bias_pct[i] - df$rmse_rel_pct[i], i,
             df$rel_bias_pct[i] + df$rmse_rel_pct[i], i,
             col = cols[i], lwd = 2)
  }
  axis_labels <- sapply(df$parameter, function(p) {
    lbl <- par_labels[[p]]
    if (is.null(lbl)) as.expression(p) else lbl   # fall back to raw name
  })
  axis(2, at = seq_len(np), labels = axis_labels, las = 2, cex.axis = 1.1)

  # Coverage annotation on right margin
  cov_labels <- ifelse(is.na(df$coverage_95), "\u2014",
                       sprintf("%.0f%%", df$coverage_95))
  mtext(cov_labels, side = 4,
        at = seq_len(np), las = 2, cex = 0.75, col = "gray30",
        line = 0.5)
  mtext("95% CI\ncoverage", side = 4, at = np + 0.8,
        las = 2, cex = 0.7, col = "gray30", line = 0.5)

  legend("bottomright",
         legend = c("Relative bias (centre)",
                    "± Relative RMSE (whiskers)"),
         pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2),
         col = "black", bty = "n", cex = 0.85)
  grid(ny = NA, col = "lightgray", lty = 1)
}

# ---- Marginal vs conditional recovery comparison plot ----
#
# Two-panel dot plot: left = conditional MLE (given states),
# right = marginal MLE (spike train only). Same x-axis scale on both
# panels to make the information cost of latent state integration
# directly visible. Parameters without coverage (marginal study does
# not compute Wald CIs) show only bias and RMSE whiskers.

plot_marginal_vs_conditional <- function(cond_summary, marg_summary,
                                         main = "Conditional vs Marginal MLE Recovery") {

  par_labels <- c(
    a_p     = expression(a[p]),
    a_s     = expression(a[s]),
    sigma_p = expression(sigma[p]),
    sigma_s = expression(sigma[s]),
    mu0     = expression(mu[0]),
    mu_0    = expression(mu[0]),     # handle both naming conventions
    rho     = expression(rho),
    c_p     = expression(c[p]),
    c_s     = expression(c[s])
  )
  cols_map <- c(
    a_p     = "#0072B2", a_s     = "#56B4E9",
    sigma_p = "#D55E00", sigma_s = "#E69F00",
    mu0     = "#2C5F2E", mu_0    = "#2C5F2E",
    rho     = "#2C5F2E",
    c_p     = "#CC79A7", c_s     = "#009E73"
  )

  # Align rows by parameter name (marginal may have fewer if some failed)
  params_cond <- cond_summary$parameter
  params_marg <- marg_summary$parameter
  params_both <- intersect(params_cond, params_marg)

  cond_df <- cond_summary[match(params_both, params_cond), ]
  marg_df <- marg_summary[match(params_both, params_marg), ]
  np      <- length(params_both)

  cols <- unname(cols_map[params_both])
  cols[is.na(cols)] <- "black"

  # Shared x range across both panels
  all_vals <- c(
    cond_df$rel_bias_pct - cond_df$rmse_rel_pct,
    cond_df$rel_bias_pct + cond_df$rmse_rel_pct,
    marg_df$rel_bias_pct - marg_df$rmse_rel_pct,
    marg_df$rel_bias_pct + marg_df$rmse_rel_pct
  )
  finite_vals <- all_vals[is.finite(all_vals)]
  xlim <- if (length(finite_vals) >= 2L)
    range(finite_vals) * 1.2
  else
    c(-20, 20)   # fallback when all replications failed or produced NA
  ylim <- c(0.5, np + 0.5)

  op <- par(mfrow = c(1, 2),
            mar   = c(4.5, 6.5, 3, 1),
            oma   = c(0, 0, 2.5, 0))

  for (panel in 1:2) {
    df    <- if (panel == 1L) cond_df else marg_df
    title <- if (panel == 1L) "Conditional MLE\n(states known)"
    else            "Marginal MLE\n(spike train only)"

    plot(df$rel_bias_pct, seq_len(np), pch = 19, cex = 1.5,
         col = cols,
         xlim = xlim, ylim = ylim,
         xlab = "Relative error (%)", ylab = "",
         yaxt = "n", main = title,
         frame.plot = FALSE)
    abline(v = 0, lty = 2, col = "gray50")
    for (i in seq_len(np)) {
      segments(df$rel_bias_pct[i] - df$rmse_rel_pct[i], i,
               df$rel_bias_pct[i] + df$rmse_rel_pct[i], i,
               col = cols[i], lwd = 2)
    }

    axis_labels <- sapply(params_both, function(p) {
      lbl <- par_labels[[p]]
      if (is.null(lbl)) as.expression(p) else lbl
    })
    axis(2, at = seq_len(np), labels = axis_labels,
         las = 2, cex.axis = 1.05)

    # Right-margin annotation: coverage (cond) or n_valid (marg)
    if (np > 0L) {
      if (panel == 1L && "coverage_95" %in% names(df)) {
        mtext(sprintf("%.0f%%", df$coverage_95), side = 4,
              at = seq_len(np), las = 2, cex = 0.72, col = "gray30", line = 0.4)
        mtext("95% CI\ncoverage", side = 4, at = np + 0.85,
              las = 2, cex = 0.65, col = "gray30", line = 0.4)
      } else {
        mtext(sprintf("n=%d", df$n_valid), side = 4,
              at = seq_len(np), las = 2, cex = 0.72, col = "gray30", line = 0.4)
        mtext("valid\nreps", side = 4, at = np + 0.85,
              las = 2, cex = 0.65, col = "gray30", line = 0.4)
      }
    }
    grid(ny = NA, col = "lightgray", lty = 1)
  }

  mtext(main, outer = TRUE, cex = 1.05, font = 2)
  par(op)
}

# ---- Profile likelihood 2x3 panel ----
#
# Plots all six profile likelihoods on a common normalized scale.
# The horizontal dashed line marks the chi-squared 95% CI threshold
# (profile LL_max - 1.92).

plot_profiles <- function(profiles,
                           true_vals = NULL,
                           main = "Profile Likelihoods") {
  par_order <- intersect(
    c("a_p", "a_s", "sigma_p", "sigma_s", "c_p", "c_s", "mu0", "rho",
      "a_ps", "a_sp"),
    names(profiles))
  nc_pl <- min(4L, length(par_order))   # 4 columns fits 8 params in 2 rows, 10 in 3
  nr_pl <- ceiling(length(par_order) / nc_pl)
  op <- par(mfrow = c(nr_pl, nc_pl), mar = c(4, 4.5, 2.5, 1), oma = c(0, 0, 2, 0))

  for (p in par_order) {
    pf      <- profiles[[p]]
    ll_norm <- pf$ll - pf$ll_mle       # normalize to 0 at MLE
    finite  <- is.finite(ll_norm)
    thresh  <- -1.92                   # chi^2(1, 0.95) / 2

    plot(pf$grid[finite], ll_norm[finite],
         type = "l", lwd = 2.5, col = "#0072B2",
         xlab = pf$label, ylab = "Profile log-lik (normalized)",
         main = "", frame.plot = FALSE,
         ylim = c(min(c(ll_norm[finite], -4), na.rm = TRUE), 0.5))

    abline(h = 0,      lty = 2, col = "gray50")
    abline(h = thresh, lty = 3, col = "#D55E00", lwd = 1.5)
    abline(v = pf$mle, lty = 2, col = "gray50")

    if (!is.null(true_vals) && p %in% names(true_vals)) {
      abline(v = true_vals[p], lty = 1, col = "#2C5F2E", lwd = 1.5)
    }

    # Approximate 95% CI from profile
    ci_idx  <- which(ll_norm >= thresh)
    n_grid  <- length(pf$grid)
    if (length(ci_idx) > 1) {
      ci_lo       <- pf$grid[min(ci_idx)]
      ci_hi       <- pf$grid[max(ci_idx)]
      lo_at_edge  <- min(ci_idx) == 1L
      hi_at_edge  <- max(ci_idx) == n_grid
      if (!lo_at_edge && !hi_at_edge) {
        rug(c(ci_lo, ci_hi), col = "#D55E00", lwd = 2, ticksize = 0.05)
      } else {
        # CI extends beyond the profiling grid: mark with open arrows instead
        if (!lo_at_edge) rug(ci_lo, col = "#D55E00", lwd = 2, ticksize = 0.05)
        if (!hi_at_edge) rug(ci_hi, col = "#D55E00", lwd = 2, ticksize = 0.05)
        mtext(sprintf("%s CI extends beyond grid",
                      if (lo_at_edge && hi_at_edge) "Both"
                      else if (lo_at_edge) "Left" else "Right"),
              side = 3, line = -1.2, cex = 0.65, col = "#D55E00", adj = 0.02)
      }
    }
    grid(col = "lightgray", lty = 1)
  }
  mtext(main, outer = TRUE, cex = 1.05, font = 2)
  par(op)
}

# ---- Marginal recovery experiment ----
# Full pipeline: simulate → pp_mle (marginal MLE from spike train) → filter
# This is the key validation of the fully marginal identifiability claim.

marginal_recovery_one <- function(true_params, duration = CANONICAL_DURATION, dt = 0.005,
                                  input_fn  = CANONICAL_INPUT_FN,
                                  estimator    = pp_mle_twostage,
                                  seed = NULL) {
  sim_res <- sim_sde_ig(duration, dt, true_params, input_fn, seed = seed)

  # c_p and c_s are fixed known constants: pass the same input_fn to the estimator
  # so the UKF prediction step removes the known input contribution.
  # estimator defaults to pp_mle; pass pp_mle_twostage for the two-stage pipeline.
  mle_result <- tryCatch(
    estimator(sim_res$spikes, true_params,
              input_fn = input_fn,
              verbose  = FALSE),
    error = function(e) NULL
  )

  # Accept code 0 (success) only. Code 1 (maxit reached without gradient criterion)
  # produces parameter vectors that may be far from the optimum; with maxit = 2000
  # and factr = 1e7 in pp_mle, code-1 exits indicate a genuine convergence failure
  # rather than a near-converged solution.
  if (is.null(mle_result) || is.na(mle_result$ll) ||
      mle_result$convergence != 0L) return(NULL)

  flt      <- mle_result$filter
  flt_grid <- filter_to_grid(flt, sim_res$time)

  delta_err  <- flt_grid$delta - sim_res$delta
  rmse_delta <- sqrt(mean(delta_err^2, na.rm = TRUE))

  fp_hat  <- mle_result$free_hat
  fp_true <- true_params$free

  # Standardize to "mu0" to match conditional MLE naming convention
  param_names <- c("a_p", "a_s", "sigma_p", "sigma_s", "mu0", "rho")
  true_vec <- c(fp_true$a_p, fp_true$a_s, fp_true$sigma_p, fp_true$sigma_s,
                fp_true$mu_0, fp_true$rho)
  hat_vec  <- c(fp_hat$a_p,  fp_hat$a_s,  fp_hat$sigma_p,  fp_hat$sigma_s,
                fp_hat$mu_0, fp_hat$rho)
  names(true_vec) <- names(hat_vec) <- param_names

  list(
    true       = true_vec,
    hat        = hat_vec,
    rmse_delta = rmse_delta,
    ll_marginal = mle_result$ll,
    n_beats    = flt$n_beats
  )
}

marginal_recovery_study <- function(N = 100L, true_params,
                                    duration = CANONICAL_DURATION,
                                    dt = 0.005,
                                    input_fn  = CANONICAL_INPUT_FN,
                                    estimator    = pp_mle_twostage,
                                    use_parallel = TRUE) {
  set.seed(RECOVERY_MASTER_SEED + 1L)   # distinct from conditional study
  seeds <- sample.int(1e6L, N)
  run_one <- function(i)
    tryCatch(marginal_recovery_one(true_params, duration, dt, input_fn,
                                   estimator = estimator,
                                   seed = seeds[i]),
             error = function(e) NULL)

  if (use_parallel && .Platform$OS.type == "unix" &&
      requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
    results <- parallel::mclapply(seq_len(N), run_one,
                                  mc.cores = n_cores, mc.set.seed = FALSE)
  } else {
    results <- lapply(seq_len(N), run_one)
  }
  results <- Filter(Negate(is.null), results)

  param_names <- c("a_p", "a_s", "sigma_p", "sigma_s", "mu0", "rho")
  summary_df <- as.data.frame(data.table::rbindlist(lapply(param_names, function(p) {
    hat_v  <- vapply(results, function(r) r$hat[[p]], numeric(1L))
    true_v <- vapply(results, function(r) r$true[[p]], numeric(1L))
    ok     <- is.finite(hat_v) & is.finite(true_v)
    if (!any(ok)) return(NULL)
    t_bar  <- mean(true_v[ok])
    h_bar  <- mean(hat_v[ok])
    bias   <- h_bar - t_bar
    rmse   <- sqrt(mean((hat_v[ok] - true_v[ok])^2))
    data.frame(parameter    = p,
               true_value   = t_bar,
               mean_hat     = h_bar,           # added
               bias         = bias,
               rel_bias_pct = 100 * bias / abs(t_bar),
               rmse         = rmse,
               rmse_rel_pct = 100 * rmse / abs(t_bar),
               n_valid      = sum(ok),
               stringsAsFactors = FALSE)
  })))

  summary_df <- as.data.frame(summary_df)

  rmse_delta_vec <- vapply(results, `[[`, numeric(1L), "rmse_delta")
  list(results    = results,
       summary    = summary_df,
       rmse_delta = rmse_delta_vec,
       N_total    = N, N_valid = length(results))
}
