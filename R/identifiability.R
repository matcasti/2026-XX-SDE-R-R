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

# ---- Single recovery experiment ----
#
# Simulates one dataset from true_params, runs the conditional MLE,
# and returns true values, estimates, and SEs side-by-side.
#
# NOTE: this is *conditional* identifiability — parameters are estimated
# given the true state trajectory. Full marginal recovery (from spike
# times alone) requires the filter in Phase 4.

single_recovery <- function(true_params,
                             duration = 300, dt = 0.005,
                             input_fn = make_double_logistic(120, 180),
                             seed = NULL) {
  sim_res <- sim_sde_ig(duration, dt, true_params, input_fn, seed = seed)
  mle     <- full_conditional_mle(sim_res)

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
                           duration = 300,
                           dt = 0.005,
                           input_fn = make_double_logistic(120, 180),
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

  summary_df <- do.call(rbind, lapply(seq_along(params), function(j) {
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
      n_valid       = n_ok,
      stringsAsFactors = FALSE
    )
  }))

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
# Displays relative bias (centre) ± relative RMSE (whiskers) per parameter.
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
  xlim <- range(all_vals, na.rm = TRUE) * 1.2
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
    ci_idx <- which(ll_norm >= thresh)
    if (length(ci_idx) > 1) {
      ci_lo <- pf$grid[min(ci_idx)]
      ci_hi <- pf$grid[max(ci_idx)]
      rug(c(ci_lo, ci_hi), col = "#D55E00", lwd = 2, ticksize = 0.05)
    }
    grid(col = "lightgray", lty = 1)
  }
  mtext(main, outer = TRUE, cex = 1.05, font = 2)
  par(op)
}

# ---- Marginal recovery experiment ----
# Full pipeline: simulate → pp_mle (marginal MLE from spike train) → filter
# This is the key validation of the fully marginal identifiability claim.

marginal_recovery_one <- function(true_params, duration = 300, dt = 0.005,
                                  input_fn = make_double_logistic(120, 180),
                                  seed = NULL) {
  sim_res <- sim_sde_ig(duration, dt, true_params, input_fn, seed = seed)

  # Inference ignores the external input: pp_mle always uses c_p = c_s = 0
  # and attributes all dynamics to the OU spectral structure.
  mle_result <- tryCatch(
    pp_mle(sim_res$spikes, true_params,
           input_fn = function(t) 0,   # no external input seen by the filter
           verbose  = FALSE),
    error = function(e) NULL
  )

  if (is.null(mle_result) || mle_result$convergence != 0L) return(NULL)

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

marginal_recovery_study <- function(N = 100L, true_params, duration = 300,
                                    dt = 0.005,
                                    input_fn = make_double_logistic(120, 180),
                                    use_parallel = TRUE) {
  set.seed(RECOVERY_MASTER_SEED + 1L)   # distinct from conditional study
  seeds <- sample.int(1e6L, N)
  run_one <- function(i)
    tryCatch(marginal_recovery_one(true_params, duration, dt, input_fn,
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
  summary_df <- do.call(rbind, lapply(param_names, function(p) {
    hat_v  <- sapply(results, function(r) r$hat[p])
    true_v <- sapply(results, function(r) r$true[p])
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
  }))

  rmse_delta_vec <- sapply(results, `[[`, "rmse_delta")
  list(results    = results,
       summary    = summary_df,
       rmse_delta = rmse_delta_vec,
       N_total    = N, N_valid = length(results))
}

# ---- Recording-duration sensitivity study ----
#
# Tests how parameter recovery (bias, RMSE) and filter RMSE on Delta(t) depend
# on recording length for either the conditional or marginal pipeline.
#
# Arguments:
#   durations    : numeric vector of recording lengths in seconds (e.g. 120 to 1800)
#   true_params  : make_model_params() output; must have c_p = c_s = 0 so that
#                  theta_true = theta* and recovery measures estimator precision only.
#   N_per_dur    : Monte Carlo replications per duration (20–50 for screening,
#                  100+ for publication-quality curves)
#   dt           : simulation time step
#   mode         : "conditional" — recovery given true states (fast, isolates
#                    estimation from latent-state integration);
#                  "marginal"    — full pipeline from spike train alone (slow but
#                    operationally relevant; requires the UKF LL fix in filter.R)
#   use_parallel : forwarded to the underlying study function
#
# Returns a list with:
#   summary_df : data.frame — one row per (duration × parameter), columns
#                  (duration, parameter, true_value, mean_hat, bias,
#                   rel_bias_pct, rmse, rmse_rel_pct, coverage_95, n_valid)
#                  For marginal mode an additional row with parameter =
#                  "delta_rmse_median" records the median filter RMSE on Delta(t).
#   N_per_dur, mode, durations, true_params : input arguments, for plotting

duration_recovery_study <- function(durations    = c(120, 300, 600, 1200),
                                    true_params,
                                    N_per_dur    = 50L,
                                    dt           = 0.005,
                                    mode         = c("conditional", "marginal"),
                                    use_parallel = TRUE) {
  mode <- match.arg(mode)

  per_dur <- lapply(durations, function(dur) {
    message(sprintf(
      "duration_recovery_study: %g s (%s mode, N = %d)",
      dur, mode, N_per_dur))

    if (mode == "conditional") {
      rec <- recovery_study(N            = N_per_dur,
                            true_params  = true_params,
                            duration     = dur,
                            dt           = dt,
                            input_fn     = function(t) 0,
                            use_parallel = use_parallel)
      df      <- rec$summary
      n_valid <- rec$N_valid

    } else {
      rec <- marginal_recovery_study(N            = N_per_dur,
                                     true_params  = true_params,
                                     duration     = dur,
                                     dt           = dt,
                                     input_fn     = function(t) 0,
                                     use_parallel = use_parallel)
      df      <- rec$summary
      n_valid <- rec$N_valid

      # Append a synthetic row for the median filter RMSE on Delta(t)
      rmse_v <- rec$rmse_delta[is.finite(rec$rmse_delta)]
      if (length(rmse_v) > 0L) {
        # Build the extra row using only the columns that exist in df,
        # which comes from marginal_recovery_study() and has no coverage_95.
        extra_row <- data.frame(
          parameter    = "delta_rmse_median",
          true_value   = NA_real_,
          mean_hat     = NA_real_,
          bias         = NA_real_,
          rel_bias_pct = NA_real_,
          rmse         = median(rmse_v),
          rmse_rel_pct = NA_real_,
          n_valid      = length(rmse_v),
          stringsAsFactors = FALSE
        )
        # Add any columns present in df but absent in extra_row as NA
        missing_cols <- setdiff(names(df), names(extra_row))
        for (col in missing_cols) extra_row[[col]] <- NA_real_
        df <- rbind(df, extra_row[, names(df)])
      }
    }

    df$duration <- dur
    df$n_valid_total <- n_valid
    df
  })

  list(
    summary_df  = do.call(rbind, per_dur),
    N_per_dur   = N_per_dur,
    mode        = mode,
    durations   = durations,
    true_params = true_params
  )
}

# ---- Plot: relative RMSE vs recording duration ----
#
# Produces a log-log line plot (one line per parameter) showing how estimation
# precision scales with recording length. Overlays a -1/2 reference slope for
# comparison against the Cramér–Rao sqrt(T) scaling expected under iid conditions.
#
# Arguments:
#   dur_study       : output of duration_recovery_study()
#   params_to_show  : parameter names to plot (subset of dur_study summary)
#   show_delta_rmse : if TRUE and mode == "marginal", overlay the Delta filter RMSE
#   main            : plot title (auto-generated if NULL)

plot_duration_study <- function(dur_study,
                                params_to_show  = c("a_p","a_s","sigma_p",
                                                    "sigma_s","mu0","rho"),
                                show_delta_rmse = TRUE,
                                main            = NULL) {
  df   <- dur_study$summary_df
  mode <- dur_study$mode
  if (is.null(main))
    main <- sprintf("Recovery vs Recording Duration — %s MLE",
                    if (mode == "conditional") "Conditional" else "Marginal")

  par_labels <- c(
    a_p              = expression(a[p]),
    a_s              = expression(a[s]),
    sigma_p          = expression(sigma[p]),
    sigma_s          = expression(sigma[s]),
    mu0              = expression(mu[0]),
    rho              = expression(rho),
    a_ps             = expression(a[ps]),
    a_sp             = expression(a[sp]),
    delta_rmse_median = expression("Filter RMSE "~hat(Delta)(t))
  )
  base_cols <- c("#0072B2","#56B4E9","#D55E00","#E69F00",
                 "#2C5F2E","#CC79A7","#009E73","#999999")

  show_pars <- intersect(params_to_show, unique(df$parameter))
  if (show_delta_rmse && "delta_rmse_median" %in% df$parameter)
    show_pars <- c(show_pars, "delta_rmse_median")

  np   <- length(show_pars)
  cols <- setNames(
    c(base_cols[seq_len(min(np, 8L))],
      rep("black", max(np - 8L, 0L))),
    show_pars
  )
  durs <- sort(unique(df$duration))

  # Collect y values for axis range
  yvals <- df$rmse_rel_pct[df$parameter %in% setdiff(show_pars, "delta_rmse_median")]
  if (show_delta_rmse && "delta_rmse_median" %in% show_pars) {
    # delta RMSE is in absolute seconds; convert to % of sd(true delta) for comparability
    # Rough reference: sd(true delta) ~ sqrt(sigma_p^2/(2*a_p) + sigma_s^2/(2*a_s))
    fp <- dur_study$true_params$free
    sd_delta <- sqrt(fp$sigma_p^2 / (2 * fp$a_p) + fp$sigma_s^2 / (2 * fp$a_s))
    df_d      <- df[df$parameter == "delta_rmse_median", ]
    delta_pct <- 100 * df_d$rmse / sd_delta
    yvals     <- c(yvals, delta_pct)
  }
  yvals <- yvals[is.finite(yvals) & yvals > 0]
  ylim  <- if (length(yvals) > 0L) range(yvals) * c(0.5, 2) else c(0.1, 200)

  par(mar = c(4.5, 5, 3, 2))
  plot(NA,
       xlim = range(durs),
       ylim = ylim,
       log  = "xy",
       xlab = "Recording duration (s)",
       ylab = "Relative RMSE (%)",
       main = main,
       frame.plot = FALSE)

  # Reference -1/2 slope (sqrt(T) CR bound)
  t_ref <- seq(min(durs), max(durs), length.out = 50)
  ref_y <- 5 * (t_ref / min(durs))^(-0.5)   # anchored at 5% for min duration
  lines(t_ref, ref_y, lty = 3, col = "gray60", lwd = 1.5)
  text(max(t_ref) * 0.85, min(ref_y) * 1.3, expression("CR slope"~T^{-1/2}),
       adj = 0, cex = 0.75, col = "gray50")

  abline(h = c(1, 5, 10), lty = 1, col = "lightgray")

  for (i in seq_along(show_pars)) {
    p <- show_pars[i]
    if (p == "delta_rmse_median") {
      if (!show_delta_rmse) next
      pd  <- df[df$parameter == p, ]
      pd  <- pd[order(pd$duration), ]
      yy  <- 100 * pd$rmse / sd_delta
      lines(pd$duration, yy, col = cols[p], lwd = 2.5,
            type = "b", pch = 17, lty = 2)
    } else {
      pd <- df[df$parameter == p, ]
      pd <- pd[order(pd$duration), ]
      ok <- is.finite(pd$rmse_rel_pct) & pd$rmse_rel_pct > 0
      if (!any(ok)) next
      lines(pd$duration[ok], pd$rmse_rel_pct[ok], col = cols[p], lwd = 2,
            type = "b", pch = 19)
    }
  }

  legend_labels <- sapply(show_pars, function(p) {
    lbl <- par_labels[[p]]
    if (is.null(lbl)) as.expression(p) else lbl
  })
  legend("topright", legend = legend_labels,
         col = unname(cols), lwd = 2,
         pch = ifelse(show_pars == "delta_rmse_median", 17, 19),
         lty = ifelse(show_pars == "delta_rmse_median", 2, 1),
         bty = "n", cex = 0.85)
  grid(col = "lightgray", lty = 1)
}
