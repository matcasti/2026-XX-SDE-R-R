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
                             input_fn = function(t) as.numeric(t >= 120 & t < 180),
                             seed = NULL) {
  sim_res <- sim_sde_ig(duration, dt, true_params, input_fn, seed = seed)
  mle     <- full_conditional_mle(sim_res)

  fp <- true_params$free
  true_vec <- c(sigma_p = fp$sigma_p, sigma_s = fp$sigma_s,
                mu0     = fp$mu_0,    kappa   = fp$kappa,
                c_p     = fp$c_p,     c_s     = fp$c_s)
  hat_vec  <- c(sigma_p = mle$sigma_p, sigma_s = mle$sigma_s,
                mu0     = mle$mu0,     kappa   = mle$kappa,
                c_p     = mle$c_p,     c_s     = mle$c_s)
  se_vec   <- c(sigma_p = mle$sigma_p_se, sigma_s = mle$sigma_s_se,
                mu0     = mle$mu0_se,     kappa   = mle$kappa_se,
                c_p     = mle$c_p_se,     c_s     = mle$c_s_se)

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

recovery_study <- function(N             = 200L,
                            true_params,
                            duration      = 300,
                            dt            = 0.005,
                            input_fn      = function(t) as.numeric(t >= 120 & t < 180),
                            use_parallel  = TRUE) {
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

  params <- c("sigma_p", "sigma_s", "mu0", "kappa", "c_p", "c_s")
  labels <- c(expression(sigma[p]), expression(sigma[s]),
              expression(mu[0]),    expression(kappa),
              expression(c[p]),     expression(c[s]))

  summary_df <- do.call(rbind, lapply(seq_along(params), function(j) {
    p      <- params[j]
    true_v <- sapply(results, function(r) r$true[p])
    hat_v  <- sapply(results, function(r) r$hat[p])
    se_v   <- sapply(results, function(r) r$se[p])
    ok     <- is.finite(hat_v) & is.finite(se_v) & se_v > 0

    n_ok   <- sum(ok)
    if (n_ok == 0L) {
      return(data.frame(
        parameter = p, true_value = NA_real_, mean_hat = NA_real_,
        bias = NA_real_, rel_bias_pct = NA_real_,
        rmse = NA_real_, rmse_rel_pct = NA_real_,
        coverage_95 = NA_real_, n_valid = 0L,
        stringsAsFactors = FALSE))
    }
    t_bar  <- mean(true_v[ok])
    h_bar  <- mean(hat_v[ok])
    bias   <- h_bar - t_bar
    rmse   <- sqrt(mean((hat_v[ok] - true_v[ok])^2))

    # 95% Wald CI coverage  (asymptotic normality of MLE)
    ci_lo  <- hat_v[ok] - 1.96 * se_v[ok]
    ci_hi  <- hat_v[ok] + 1.96 * se_v[ok]
    cover  <- mean(ci_lo <= true_v[ok] & true_v[ok] <= ci_hi)

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
    sigma_p = expression(sigma[p]),
    sigma_s = expression(sigma[s]),
    mu0     = expression(mu[0]),
    kappa   = expression(kappa),
    c_p     = expression(c[p]),
    c_s     = expression(c[s])
  )

  np   <- nrow(df)
  cols <- c("#0072B2", "#0072B2", "#D55E00", "#D55E00", "#2C5F2E", "#2C5F2E")

  par(mar = c(4.5, 6, 3, 5))
  plot(df$rel_bias_pct, seq_len(np), pch = 19, cex = 1.5,
       col  = cols,
       xlim = range(c(df$rel_bias_pct - df$rmse_rel_pct,
                      df$rel_bias_pct + df$rmse_rel_pct)) * 1.15,
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
  mtext(sprintf("%.0f%%", df$coverage_95), side = 4,
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

# ---- Profile likelihood 2x3 panel ----
#
# Plots all six profile likelihoods on a common normalized scale.
# The horizontal dashed line marks the chi-squared 95% CI threshold
# (profile LL_max - 1.92).

plot_profiles <- function(profiles,
                           true_vals = NULL,
                           main = "Profile Likelihoods") {
  par_order <- c("mu0", "kappa", "sigma_p", "sigma_s", "c_p", "c_s")
  op <- par(mfrow = c(2, 3), mar = c(4, 4.5, 2.5, 1), oma = c(0, 0, 2, 0))

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
