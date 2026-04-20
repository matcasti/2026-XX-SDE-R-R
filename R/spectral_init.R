# ============================================================
# R/spectral_init.R
#
# Data-driven starting-value estimation for pp_mle().
# Exploits the fact that standard HRV statistics have closed-form
# expressions in terms of the SDE-IG free parameters, giving a
# cheap, analytically grounded initialization that is tailored to
# the actual recording rather than to literature defaults.
#
# Core identities used (all derived from the SDE-IG generative model):
#
#   E[tau]       = mu_0 * exp(sigma_d2 / 2)          [log-normal offset]
#   Corr[tau_j]  ≈ wp*exp(-a_p*j*mu_0) + ws*exp(-a_s*j*mu_0)
#                                                    [biexponential ACF]
#   sigma_d2     = sigma_p^2/(2*a_p) + sigma_s^2/(2*a_s)
#   kappa        = N / sum(1/tau_k - 1/mu_obs)        [stationary IG MLE]
#
# The ACF biexponential fit recovers (a_p, a_s, wp, ws);
# the spectral weights give (sigma_p, sigma_s);
# the mean gives mu_0;
# the harmonic residual gives kappa and rho.
#
# Returns a make_model_params()-compatible object suitable as
# the params_init argument of pp_mle().
#
# Requires: model_functions.R sourced first.
# ============================================================

# ---- Biexponential ACF fitting ----
# Internal helper: fit r(j) = wp*exp(-dp*j) + (1-wp)*exp(-ds*j)
# where dp = a_p*mu_obs and ds = a_s*mu_obs are dimensionless.
# Returns list(wp, ws, dp, ds) or NULL on failure.

.fit_biexp_acf <- function(acf_v, max_lag) {
  lags <- seq_len(max_lag)

  # Primary attempt: port-constrained NLS
  fit <- tryCatch(
    nls(
      acf_v ~ wp * exp(-dp * lags) + (1 - wp) * exp(-ds * lags),
      start     = list(wp = 0.25, dp = 1.50, ds = 0.15),
      lower     = c(wp = 0.01, dp = 0.10, ds = 0.005),
      upper     = c(wp = 0.99, dp = 20.0,  ds = 2.00),
      algorithm = "port",
      control   = nls.control(maxiter = 300L, tol = 1e-6, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )

  # Fallback: Nelder-Mead over the same parameterisation
  if (is.null(fit)) {
    obj <- function(th) {
      wp <- plogis(th[1L]); dp <- exp(th[2L]); ds <- exp(th[3L])
      sum((acf_v - wp * exp(-dp * lags) - (1 - wp) * exp(-ds * lags))^2)
    }
    opt <- tryCatch(
      optim(c(qlogis(0.25), log(1.5), log(0.15)), obj,
            method = "Nelder-Mead",
            control = list(maxit = 2000L, reltol = 1e-8)),
      error = function(e) list(value = Inf)
    )
    if (is.finite(opt$value)) {
      cf <- opt$par
      return(list(wp = plogis(cf[1L]), ws = 1 - plogis(cf[1L]),
                  dp = exp(cf[2L]),     ds = exp(cf[3L]),
                  source = "nelder_mead"))
    }
    return(NULL)
  }

  cf <- coef(fit)
  list(wp = cf[["wp"]], ws = 1 - cf[["wp"]],
       dp = cf[["dp"]], ds = cf[["ds"]],
       source = "nls")
}

# ---- Main spectral initialisation function ----

spectral_init <- function(rr_vec, max_lag = 20L, verbose = FALSE) {
  stopifnot(
    is.numeric(rr_vec), length(rr_vec) >= max_lag + 5L,
    all(rr_vec > 0)
  )

  mu_obs <- mean(rr_vec)
  sdnn2  <- var(rr_vec)

  # ----------------------------------------------------------------
  # Step 1 — biexponential ACF fit
  # ACF lags 1:max_lag are approximately:
  #   r(j) = wp*exp(-a_p*mu_obs*j) + ws*exp(-a_s*mu_obs*j)
  # The per-lag decay rates dp = a_p*mu_obs, ds = a_s*mu_obs are
  # dimensionless; dividing by mu_obs recovers the Hz rates.
  # ----------------------------------------------------------------
  acf_v <- as.numeric(acf(rr_vec, lag.max = max_lag, plot = FALSE)$acf)[-1L]
  biexp <- .fit_biexp_acf(acf_v, max_lag)

  if (!is.null(biexp)) {
    wp_hat <- biexp$wp; ws_hat <- biexp$ws
    a_p_hat <- biexp$dp / mu_obs
    a_s_hat <- biexp$ds / mu_obs
    # Enforce a_p > a_s: vagal branch has faster decay
    if (a_p_hat < a_s_hat) {
      a_p_hat <- biexp$ds / mu_obs
      a_s_hat <- biexp$dp / mu_obs
      wp_hat  <- biexp$ws
      ws_hat  <- biexp$wp
    }
    if (verbose) message(sprintf(
      "spectral_init: ACF fit (%s) a_p=%.3f a_s=%.3f wp=%.3f",
      biexp$source, a_p_hat, a_s_hat, wp_hat))
  } else {
    # Literature fallback — identical to make_model_params() defaults
    a_p_hat <- 2.0; a_s_hat <- 0.2
    wp_hat  <- 0.25; ws_hat  <- 0.75
    if (verbose) message("spectral_init: ACF fit failed; using literature defaults")
  }

  # Clamp to physiologically plausible range
  a_p_hat <- max(min(a_p_hat, 15.0), 0.5)
  a_s_hat <- max(min(a_s_hat, 1.0),  0.02)
  if (a_p_hat <= a_s_hat) a_p_hat <- a_s_hat * 5

  # ----------------------------------------------------------------
  # Step 2 — spectral variances from SDNN and ACF weights
  # Var[tau] ≈ mu_0^2 * sigma_d2  (OU contribution; IG intrinsic is
  # secondary for typical rho).  Hence sigma_d2 ≈ SDNN^2 / mu_obs^2.
  # The ACF weights distribute the total variance between branches:
  #   sigma_p^2/(2*a_p) = wp * sigma_d2
  #   sigma_s^2/(2*a_s) = ws * sigma_d2
  # ----------------------------------------------------------------
  sigma_d2    <- sdnn2 / (mu_obs^2)      # dimensionless
  sigma_p_hat <- sqrt(2 * a_p_hat * wp_hat * sigma_d2)
  sigma_s_hat <- sqrt(2 * a_s_hat * ws_hat * sigma_d2)
  sigma_p_hat <- max(sigma_p_hat, 0.01)
  sigma_s_hat <- max(sigma_s_hat, 0.01)

  # ----------------------------------------------------------------
  # Step 3 — mu_0 with log-normal mean correction
  # E[tau] = mu_0 * exp(sigma_d2 / 2), so
  # mu_0 = mu_obs * exp(-sigma_d2 / 2)
  # ----------------------------------------------------------------
  mu_0_hat <- mu_obs * exp(-sigma_d2 / 2)
  mu_0_hat <- max(mu_0_hat, 0.30)

  # ----------------------------------------------------------------
  # Step 4 — rho from the stationary IG MLE
  # Treating Delta_k ≈ 0 (ignoring OU variation for this step):
  #   kappa_MLE = 1 / (mean(1/tau) - 1/mean(tau))
  # This estimator is biased upward when sigma_d2 > 0, but the bias
  # is absorbed by the subsequent pp_mle optimisation.
  # ----------------------------------------------------------------
  harm_excess <- mean(1 / rr_vec) - 1 / mu_obs   # >= 0 by Jensen
  kappa_hat   <- if (harm_excess > 1e-8) 1 / harm_excess else 10.0
  kappa_hat   <- max(min(kappa_hat, 1e4), 0.5)
  rho_hat     <- sqrt(mu_0_hat / kappa_hat)
  rho_hat     <- max(min(rho_hat, 0.80), 0.05)

  if (verbose) message(sprintf(
    "spectral_init: mu_0=%.3f rho=%.3f sigma_p=%.3f sigma_s=%.3f",
    mu_0_hat, rho_hat, sigma_p_hat, sigma_s_hat))

  make_model_params(
    a_p     = a_p_hat,   a_s     = a_s_hat,
    sigma_p = sigma_p_hat, sigma_s = sigma_s_hat,
    mu_0    = mu_0_hat,  rho     = rho_hat
  )
}

# ---- Predicted moments from model parameters ----
# Returns the six HRV moments that are analytically constrained by
# the SDE-IG free parameters.  Used by check_moment_consistency()
# in mle.R.

predicted_hrv_moments <- function(a_p, a_s, sigma_p, sigma_s, mu_0, rho) {
  kappa    <- kappa_from_rho(mu_0, rho)
  sigma_d2 <- sigma_p^2 / (2 * a_p) + sigma_s^2 / (2 * a_s)

  # E[tau]
  e_tau   <- mu_0 * exp(sigma_d2 / 2)

  # Var[tau] = IG intrinsic + OU-driven
  ig_var  <- mu_0^3 * exp(9 * sigma_d2 / 2) / kappa
  ou_var  <- mu_0^2 * exp(sigma_d2) * expm1(sigma_d2)
  var_tau <- ig_var + ou_var

  # RMSSD^2 = 2 * Var[tau] * (1 - r(1))
  # r(1) ≈ wp*exp(-a_p*mu_0) + ws*exp(-a_s*mu_0)
  wp      <- (sigma_p^2 / (2 * a_p)) / sigma_d2
  ws      <- 1 - wp
  r1      <- wp * exp(-a_p * e_tau) + ws * exp(-a_s * e_tau)
  rmssd2  <- 2 * var_tau * (1 - r1)

  # LF and HF band powers: integral of Lorentzian S_p(f) + S_s(f)
  # over the standard HRV bands
  lorentz_integral <- function(a, f_lo, f_hi) {
    (1 / pi) * (atan(2 * pi * f_hi / a) - atan(2 * pi * f_lo / a))
  }
  hf_power <- sigma_p^2 / (2 * a_p) * lorentz_integral(a_p, 0.15, 0.40) +
    sigma_s^2 / (2 * a_s) * lorentz_integral(a_s, 0.15, 0.40)
  lf_power <- sigma_p^2 / (2 * a_p) * lorentz_integral(a_p, 0.04, 0.15) +
    sigma_s^2 / (2 * a_s) * lorentz_integral(a_s, 0.04, 0.15)

  list(
    e_tau    = e_tau,
    sdnn     = sqrt(var_tau),
    rmssd    = sqrt(max(rmssd2, 0)),
    lf_power = lf_power,
    hf_power = hf_power,
    lf_hf    = lf_power / max(hf_power, 1e-12)
  )
}
