# ============================================================
# R/sensitivity.R
#
# A3.3 — Decay-rate sensitivity analysis.
# Tests whether conclusions about identifiability and filter
# performance are stable across physiologically plausible
# variation in the free decay-rate parameters a_p and a_s.
#
# Design: a_p and a_s are each varied independently over a
# ±50% grid around their reference values. For each grid point
# we (1) re-run the conditional MLE recovery study (N = 50
# replications, shorter than the main study for speed) and
# record bias and RMSE for sigma_p, sigma_s, mu_0, kappa;
# (2) re-run the filter on the reference simulation and record
# filter RMSE on Delta(t) and total log-likelihood.
#
# Requires: model_functions.R, mle.R, identifiability.R,
#           filter.R sourced first.
# ============================================================

SENSITIVITY_SEED <- 9999L

# ---- Grid of structural parameter values ----
# Reference: a_p = 2.0 Hz, a_s = 0.2 Hz.
# Range: ±50% in multiplicative steps.

sensitivity_grid <- function(ref_ap = 2.0,
                             ref_as = 0.2,
                             multipliers = c(0.5, 0.75, 1.0, 1.25, 1.5)) {
  expand.grid(
    ap_mult = multipliers,
    as_mult = multipliers
  ) |>
    transform(
      a_p = ref_ap * ap_mult,
      a_s = ref_as * as_mult
    ) |>
    subset(a_p > a_s)   # structural constraint: vagal must be faster
}

# ---- Single grid-point evaluation ----
# Returns a one-row data.frame of summary statistics.

sensitivity_one <- function(a_p_val, a_s_val,
                            base_params,
                            sim_res_ref,    # reference simulation (fixed seed)
                            N_rep   = 50L,
                            duration = 300,
                            dt       = 0.005,
                            input_fn = make_double_logistic(120, 180)) {

  # Build modified parameter object — vary true decay rates across the grid
  params_mod <- base_params
  params_mod$free$a_p <- a_p_val
  params_mod$free$a_s <- a_s_val

  # (1) Conditional MLE recovery — N_rep short replications
  set.seed(SENSITIVITY_SEED)
  seeds <- sample.int(1e6L, N_rep)

  results <- lapply(seq_len(N_rep), function(i) {
    tryCatch(
      single_recovery(params_mod, duration, dt, input_fn, seed = seeds[i]),
      error = function(e) NULL
    )
  })
  results <- Filter(Negate(is.null), results)

  # Aggregate bias and RMSE for the four IG + OU-amplitude parameters
  agg <- lapply(c("a_p", "a_s", "sigma_p", "sigma_s", "mu0", "rho"), function(p) {
    hat_v  <- sapply(results, function(r) r$hat[p])
    true_v <- sapply(results, function(r) r$true[p])
    ok     <- is.finite(hat_v) & is.finite(true_v)
    if (sum(ok) == 0L)
      return(c(bias = NA_real_, rmse = NA_real_, rel_rmse_pct = NA_real_))
    bias     <- mean(hat_v[ok] - true_v[ok])
    rmse     <- sqrt(mean((hat_v[ok] - true_v[ok])^2))
    t_bar    <- mean(true_v[ok])
    rel_rmse <- if (abs(t_bar) > 1e-10) 100 * rmse / abs(t_bar) else NA_real_
    c(bias = bias, rmse = rmse, rel_rmse_pct = rel_rmse)
  })
  names(agg) <- c("a_p", "a_s", "sigma_p", "sigma_s", "mu0", "rho")

  # (2) Filter RMSE on the reference simulation with a_p / a_s FIXED to grid values.
  # This is intentionally NOT a joint estimate: the goal is to isolate how much
  # filter accuracy degrades when the decay rates are *misspecified*, independently
  # of estimation variance.
  flt <- tryCatch({
    params_for_filter        <- base_params
    params_for_filter$free$a_p <- a_p_val
    params_for_filter$free$a_s <- a_s_val
    pp_ukf(sim_res_ref$spikes, params_for_filter,
           input_fn = sim_res_ref$input_fn)
  }, error = function(e) NULL)

  flt_rmse <- NA_real_
  flt_ll   <- NA_real_
  if (!is.null(flt)) {
    flt_grid <- filter_to_grid(flt, sim_res_ref$time)
    delta_err <- flt_grid$delta - sim_res_ref$delta
    flt_rmse  <- sqrt(mean(delta_err^2, na.rm = TRUE))
    flt_ll    <- flt$ll
  }

  data.frame(
    a_p               = a_p_val,
    a_s               = a_s_val,
    rel_rmse_a_p      = agg$a_p["rel_rmse_pct"],
    rel_rmse_a_s      = agg$a_s["rel_rmse_pct"],
    rel_rmse_sigma_p  = agg$sigma_p["rel_rmse_pct"],
    rel_rmse_sigma_s  = agg$sigma_s["rel_rmse_pct"],
    rel_rmse_mu0      = agg$mu0["rel_rmse_pct"],
    rel_rmse_rho      = agg$rho["rel_rmse_pct"],
    filter_rmse_delta = flt_rmse,
    filter_ll         = flt_ll,
    n_valid           = length(results),
    stringsAsFactors  = FALSE
  )
}

# ---- Full sensitivity study ----

run_sensitivity <- function(base_params,
                            sim_res_ref,
                            multipliers   = c(0.5, 0.75, 1.0, 1.25, 1.5),
                            N_rep         = 50L,
                            use_parallel  = TRUE) {
  grid <- sensitivity_grid(
    ref_ap       = base_params$free$a_p,
    ref_as       = base_params$free$a_s,
    multipliers  = multipliers
  )

  run_row <- function(i) {
    sensitivity_one(grid$a_p[i], grid$a_s[i],
                    base_params, sim_res_ref,
                    N_rep = N_rep)
  }

  if (use_parallel && .Platform$OS.type == "unix" &&
      requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
    rows    <- parallel::mclapply(seq_len(nrow(grid)), run_row,
                                  mc.cores = n_cores, mc.set.seed = FALSE)
  } else {
    rows <- lapply(seq_len(nrow(grid)), run_row)
  }

  do.call(rbind, Filter(Negate(is.null), rows))
}

# ---- Heatmap plot ----
# One panel per summary statistic, x = a_p, y = a_s.
# The reference cell (multiplier = 1.0, 1.0) is outlined in black.

plot_sensitivity <- function(sens_df,
                             ref_ap = 2.0,
                             ref_as = 0.2) {
  metrics <- c("rel_rmse_a_p",    "rel_rmse_a_s",
               "rel_rmse_sigma_p","rel_rmse_sigma_s",
               "rel_rmse_mu0",    "rel_rmse_rho",
               "filter_rmse_delta")
  labels  <- c(expression("Rel. RMSE" ~ a[p] ~ "(%)"),
               expression("Rel. RMSE" ~ a[s] ~ "(%)"),
               expression("Rel. RMSE" ~ sigma[p] ~ "(%)"),
               expression("Rel. RMSE" ~ sigma[s] ~ "(%)"),
               expression("Rel. RMSE" ~ mu[0] ~ "(%)"),
               expression("Rel. RMSE" ~ rho ~ "(%)"),
               expression("Filter RMSE" ~ hat(Delta)(t)))

  ap_vals <- sort(unique(sens_df$a_p))
  as_vals <- sort(unique(sens_df$a_s))
  n_ap    <- length(ap_vals)
  n_as    <- length(as_vals)

  op <- par(mfrow = c(2, 4), mar = c(4, 4.5, 2.5, 3),
            oma  = c(0, 0, 2.5, 0))

  for (mi in seq_along(metrics)) {
    m     <- metrics[mi]
    vals  <- matrix(NA_real_, n_as, n_ap)
    for (i in seq_len(nrow(sens_df))) {
      ci <- which(ap_vals == sens_df$a_p[i])
      ri <- which(as_vals == sens_df$a_s[i])
      if (length(ci) == 1 && length(ri) == 1)
        vals[ri, ci] <- sens_df[[m]][i]
    }
    all_na <- all(is.na(vals))
    col_range <- if (all_na) c(0, 1) else range(vals, na.rm = TRUE)
    if (all_na) {
      plot.new()
      title(main = labels[[mi]])
      text(0.5, 0.5, "no data", col = "gray50", cex = 1.2)
      next
    }
    image(ap_vals, as_vals, t(vals),
          col  = hcl.colors(20, "YlOrRd", rev = TRUE),
          zlim = col_range,
          xlab = expression(a[p] ~ "(Hz)"),
          ylab = expression(a[s] ~ "(Hz)"),
          main = labels[[mi]], axes = FALSE)
    axis(1, at = ap_vals, labels = round(ap_vals, 2))
    axis(2, at = as_vals, labels = round(as_vals, 3))
    box()
    # Outline reference cell
    ref_ci <- which.min(abs(ap_vals - ref_ap))
    ref_ri <- which.min(abs(as_vals - ref_as))
    dap    <- diff(ap_vals)[1] / 2
    das    <- diff(as_vals)[1] / 2
    rect(ap_vals[ref_ci] - dap, as_vals[ref_ri] - das,
         ap_vals[ref_ci] + dap, as_vals[ref_ri] + das,
         border = "black", lwd = 2.5)
  }
  mtext("Decay-Rate Sensitivity (Free Parameters)", outer = TRUE,
        cex = 1.05, font = 2)
  par(op)
}
