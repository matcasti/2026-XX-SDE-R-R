# ============================================================
# R/validation.R
#
# A3.1 — Real data validation on publicly available benchmarks.
# A3.2 — Benchmark comparison against Barbieri 2005 (IG-PP).
#
# Datasets (download separately; not redistributed):
#   PhysioNet Fantasia:  https://physionet.org/content/fantasia/1.0.0/
#   PhysioNet PTBDB tilt-table: https://physionet.org/content/tilt/1.0/
#
# Expected local paths (set via options or environment variables):
#   options(sde_ig.fantasia_path = "/path/to/fantasia/")
#   options(sde_ig.tilt_path    = "/path/to/tilt/")
#
# Requires: model_functions.R, mle.R, filter.R sourced first.
# ============================================================


# ---- Data loading helpers ----
# TODO (A3.1): implement WFDB RR-interval reader.
# The rhrv or RHRV package reads .qrs annotation files directly.
# Alternatively, export RR intervals to plain text via wfdb2mat and read here.

load_fantasia_rr <- function(subject_id, path = getOption("sde_ig.fantasia_path")) {
  # Returns a numeric vector of RR intervals (seconds) for one Fantasia subject.
  # Subject IDs: "f1o01" … "f1o20" (older) / "f1y01" … "f1y20" (younger).
  stop("TODO: implement WFDB reader — see RHRV::LoadBeatRHRV or wfdb toolbox.")
}

load_tilt_rr <- function(subject_id, path = getOption("sde_ig.tilt_path")) {
  # Returns list(rr = <numeric>, protocol = <data.frame with phase timestamps>)
  stop("TODO: implement tilt-table RR loader.")
}


# ---- SDE-IG fit to a real RR vector ----

fit_sde_ig <- function(rr_vec,
                       params_init = make_model_params(),
                       input_fn    = function(t) 0,
                       verbose     = FALSE) {
  # Reconstruct spike times from RR intervals
  spikes <- c(0, cumsum(rr_vec))
  pp_mle(spikes, params_init, input_fn,
         optimize_gains = FALSE, verbose = verbose)
}


# ---- Barbieri 2005 baseline: stationary IG-PP ----
# The Barbieri 2005 model is a local-maximum-likelihood point process
# with IG inter-event distribution and an AR history term.
# Minimal version for benchmarking: stationary IG renewal process
# (no history, no time variation) — the simplest nested competitor.
# Full local-likelihood extension is a TODO for a fairer comparison.

barbieri_ig_fit <- function(rr_vec) {
  # MLE for stationary IG(mu, kappa) — closed form.
  mu_hat    <- mean(rr_vec)
  kappa_hat <- length(rr_vec) / sum(1 / rr_vec - 1 / mu_hat)
  kappa_hat <- max(kappa_hat, 1e-4)

  ll <- sum(log_ig_pdf(rr_vec, mu_hat, kappa_hat))

  n_params <- 2L
  aic <- -2 * ll + 2 * n_params
  bic <- -2 * ll + log(length(rr_vec)) * n_params

  # KS statistic against IG(mu_hat, kappa_hat) CDF
  ks <- ks.test(rr_vec, function(x) {
    pnorm(sqrt(kappa_hat / x) * (x / mu_hat - 1)) +
      exp(2 * kappa_hat / mu_hat) *
      pnorm(-sqrt(kappa_hat / x) * (x / mu_hat + 1))
  })

  list(mu = mu_hat, kappa = kappa_hat,
       ll = ll, aic = aic, bic = bic,
       ks_stat = ks$statistic, ks_p = ks$p.value,
       n_params = n_params, n_beats = length(rr_vec))
}


# ---- SDE-IG model selection statistics ----

sde_ig_model_stats <- function(fit_result, rr_vec,
                               grid_dt = 0.005) {
  # fit_result: output of fit_sde_ig() (which wraps pp_mle)
  ll      <- fit_result$ll
  n_beats <- length(rr_vec)
  n_params <- 4L
  aic <- -2 * ll + 2 * n_params
  bic <- -2 * ll + log(n_beats) * n_params

  # Reconstruct a dense time grid to run compute_time_rescaling.
  # The filter only produces estimates at beat times; interpolate to a fine grid.
  spikes   <- c(0, cumsum(rr_vec))
  T_end    <- tail(spikes, 1L)
  tg       <- seq(0, T_end, by = grid_dt)
  flt_grid <- filter_to_grid(fit_result$filter, tg)

  # Build a sim_res-like list that compute_time_rescaling can consume.
  res_proxy <- list(
    time   = tg,
    spikes = spikes,
    mu     = flt_grid$mu,
    params = fit_result$params_hat
  )

  Lk  <- compute_time_rescaling(res_proxy)
  uk  <- 1 - exp(-Lk[is.finite(Lk)])
  ks  <- ks.test(uk, "punif")

  list(ll = ll, aic = aic, bic = bic,
       ks_stat = ks$statistic, ks_p = ks$p.value,
       n_params = n_params, n_beats = n_beats)
}


# ---- Head-to-head comparison table ----

compare_models <- function(sde_ig_stats, barbieri_stats) {
  data.frame(
    Model    = c("SDE-IG (this work)", "Barbieri 2005 (stationary IG)"),
    N_params = c(sde_ig_stats$n_params, barbieri_stats$n_params),
    LogLik   = sprintf("%.1f", c(sde_ig_stats$ll,    barbieri_stats$ll)),
    AIC      = sprintf("%.1f", c(sde_ig_stats$aic,   barbieri_stats$aic)),
    BIC      = sprintf("%.1f", c(sde_ig_stats$bic,   barbieri_stats$bic)),
    KS_stat  = sprintf("%.4f", c(sde_ig_stats$ks_stat, barbieri_stats$ks_stat)),
    KS_p     = sprintf("%.4f", c(sde_ig_stats$ks_p,    barbieri_stats$ks_p)),
    stringsAsFactors = FALSE
  )
}
