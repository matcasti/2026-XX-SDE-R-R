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

# ---- Power-weighted spectral band centroid ----
# Internal helper used by band_centroid_init() and spectral_init().
# Returns the power-weighted centroid frequency within [f_lo, f_hi] from
# a one-sided PSD defined on a uniform frequency grid.

.band_centroid <- function(freqs, psd, f_lo, f_hi) {
  idx <- freqs >= f_lo & freqs <= f_hi & freqs > 0
  if (!any(idx)) return((f_lo + f_hi) / 2)
  num <- sum(freqs[idx] * psd[idx])
  den <- sum(psd[idx])
  if (den < .Machine$double.eps) return((f_lo + f_hi) / 2)
  num / den
}

# ---- Band-centroid decay-rate initializer ----
# Estimates (a_p, a_s, sigma_p, sigma_s, mu_0, rho) from the tachogram
# spectrum without any nonlinear fitting. Directly analogous to
# estimateParams() in the JS client-side implementation.
#
# For an OU process with spectral density S(f) ∝ 1/(a² + (2πf)²),
# the power-weighted centroid of the HF (LF) band estimates a/(2π),
# giving a = 2π·f_centroid. More robust than biexponential ACF fitting
# for short (<300 s) or spectrally flat recordings.
#
# sigma decomposition follows from S_Δ = σ_p²/(2a_p) + σ_s²/(2a_s)
# and the empirical spectral fractions hfFrac, lfFrac:
#   σ_p² = hfFrac · σ_Δ² · 2a_p  (spectral variance allocation)
#   σ_s² = lfFrac · σ_Δ² · 2a_s

band_centroid_init <- function(rr_vec, fs = 4.0) {
  stopifnot(is.numeric(rr_vec), length(rr_vec) >= 20L, all(rr_vec > 0))

  mu_obs   <- mean(rr_vec)
  rho_obs  <- sd(rr_vec) / mu_obs     # observed CV (upper bound on σ_Δ)

  # Resample to uniform 4-Hz grid via linear interpolation
  cum_t    <- cumsum(c(0, rr_vec))
  t_rs     <- seq(0, tail(cum_t, 1L), by = 1 / fs)
  rr_rs    <- approx(cum_t[-length(cum_t)], rr_vec, xout = t_rs, rule = 2L)$y
  n_rs     <- length(rr_rs)

  # Hanning window + zero-pad to next power of 2
  w     <- 0.5 * (1 - cos(2 * pi * seq_len(n_rs) / (n_rs - 1)))
  N_fft <- nextn(n_rs * 2L, factors = 2L)
  padded <- c((rr_rs - mean(rr_rs)) * w, rep(0, N_fft - n_rs))

  sp_full <- Mod(fft(padded))[seq_len(N_fft / 2L + 1L)]^2 / (n_rs * fs)
  freqs   <- seq(0, fs / 2, length.out = N_fft / 2L + 1L)
  df      <- freqs[2L]
  # One-sided PSD scaling
  sp_full[-c(1L, length(sp_full))] <- sp_full[-c(1L, length(sp_full))] * 2

  # Decay rates from band centroids: a = 2π·f_centroid
  a_p_hat <- max(min(2 * pi * .band_centroid(freqs, sp_full, 0.15, 0.40), 5.0), 0.5)
  a_s_hat <- max(min(2 * pi * .band_centroid(freqs, sp_full, 0.04, 0.15), 1.5), 0.1)
  if (a_p_hat <= a_s_hat) a_p_hat <- a_s_hat * 5

  # LF/HF band powers → spectral fractions
  hf_pow  <- sum(sp_full[freqs >= 0.15 & freqs <= 0.40]) * df
  lf_pow  <- sum(sp_full[freqs >= 0.04 & freqs <= 0.15]) * df
  tot_pow <- max(hf_pow + lf_pow, .Machine$double.eps)
  hf_frac <- hf_pow / tot_pow
  lf_frac <- lf_pow / tot_pow

  # σ_Δ² ≈ ρ_obs² (tachogram CV is a good proxy for the OU variance)
  var_d       <- rho_obs^2
  sigma_p_hat <- sqrt(max(hf_frac * var_d * 2 * a_p_hat, 1e-6))
  sigma_s_hat <- sqrt(max(lf_frac * var_d * 2 * a_s_hat, 1e-6))

  # μ₀ with log-normal mean correction
  mu_0_hat <- max(mu_obs * exp(-var_d / 2), 0.30)

  # ρ from stationary IG MLE
  harm_excess <- mean(1 / rr_vec) - 1 / mu_obs
  kappa_hat   <- max(min(if (harm_excess > 1e-8) 1 / harm_excess else 10.0, 1e4), 0.5)
  rho_hat     <- max(min(sqrt(mu_0_hat / kappa_hat), 0.80), 0.05)

  make_model_params(
    a_p     = a_p_hat,   a_s     = a_s_hat,
    sigma_p = sigma_p_hat, sigma_s = sigma_s_hat,
    mu_0    = mu_0_hat,  rho     = rho_hat
  )
}

# ---- Bivariate AR(1) coupling initializer on band-filtered proxies ----
# Band-pass filtered RR intervals serve as surrogates for p(t) (HF) and
# s(t) (LF). A bivariate AR(1) on these proxies estimates the off-diagonal
# elements alpha_ps and alpha_sp; converting to continuous time via
# a_ps ≈ alpha_ps / Δt gives deterministic O(n) starting values for the
# coupled OU MLE. Stability is enforced by the same guard used in estimateCoupling()
# of the JS client-side code (Castillo-Aguilar 2026).

band_filtered_coupling_init <- function(rr_vec, a_p, a_s, fs = 4.0) {
  stopifnot(is.numeric(rr_vec), length(rr_vec) >= 40L)
  dt_filt <- 1 / fs

  # Resample
  cum_t <- cumsum(c(0, rr_vec))
  t_rs  <- seq(0, tail(cum_t, 1L), by = dt_filt)
  rr_rs <- approx(cum_t[-length(cum_t)], rr_vec, xout = t_rs, rule = 2L)$y
  nr    <- length(rr_rs)
  if (nr < 40L) return(list(a_ps = 0, a_sp = 0))
  x <- rr_rs - mean(rr_rs)

  # FFT band-pass filter
  bp_fft <- function(xv, f_lo, f_hi) {
    N_f    <- nextn(length(xv) * 2L, factors = 2L)
    sp_f   <- fft(c(xv, rep(0, N_f - length(xv))))
    f_vec  <- seq(0, fs, length.out = N_f)
    keep   <- (f_vec >= f_lo & f_vec <= f_hi) |
      (f_vec >= (fs - f_hi) & f_vec <= (fs - f_lo))
    sp_f[!keep] <- 0 + 0i
    Re(fft(sp_f, inverse = TRUE) / N_f)[seq_len(length(xv))]
  }

  pP <- bp_fft(x, 0.15, 0.40)   # HF proxy ≡ p(t)
  sP <- bp_fft(x, 0.04, 0.15)   # LF proxy ≡ s(t)

  m    <- nr - 1L
  pp   <- pP[seq_len(m)]; spv <- sP[seq_len(m)]
  pn   <- pP[-1L];        sn  <- sP[-1L]

  sp2   <- sum(pp^2); ss2   <- sum(spv^2); sps_c <- sum(pp * spv)
  spnp  <- sum(pp * pn); spns <- sum(spv * pn)
  ssnp  <- sum(pp * sn); ssns <- sum(spv * sn)
  det_c <- sp2 * ss2 - sps_c^2 + .Machine$double.eps

  # OLS off-diagonals of the bivariate AR(1)
  alpha_ps <- (sp2 * spns - sps_c * spnp) / det_c   # s→p  in discrete time
  alpha_sp <- (ss2 * ssnp - sps_c * ssns) / det_c   # p→s  in discrete time

  # First-order continuous-time conversion
  a_ps_hat <- alpha_ps / dt_filt
  a_sp_hat <- alpha_sp / dt_filt

  # Stability guard: |a_ps · a_sp| < 0.8 · a_p · a_s
  max_cpl <- 0.8 * a_p * a_s
  mag     <- sqrt(abs(a_ps_hat * a_sp_hat)) + .Machine$double.eps
  if (mag^2 > max_cpl) {
    sc       <- sqrt(max_cpl) / mag
    a_ps_hat <- a_ps_hat * sc;  a_sp_hat <- a_sp_hat * sc
  }
  a_ps_hat <- max(min(a_ps_hat, a_p * 0.85), 0)
  a_sp_hat <- max(min(a_sp_hat, a_s * 0.85), 0)

  list(a_ps = a_ps_hat, a_sp = a_sp_hat)
}

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

  # Band-centroid estimates: always computed as a reliability check.
  bc    <- band_centroid_init(rr_vec)
  bc_ap <- bc$free$a_p;  bc_as <- bc$free$a_s

  if (!is.null(biexp)) {
    a_p_hat <- biexp$dp / mu_obs
    a_s_hat <- biexp$ds / mu_obs
    wp_hat  <- biexp$wp;  ws_hat <- biexp$ws
    if (a_p_hat < a_s_hat) {
      a_p_hat <- biexp$ds / mu_obs;  a_s_hat <- biexp$dp / mu_obs
      wp_hat  <- biexp$ws;           ws_hat  <- biexp$wp
    }
    # Ensemble: average log-scale estimates from ACF and band-centroid.
    # This damps the ACF estimator's tendency to lock onto a single dominant
    # ACF lag when the spectrum is non-peaked (a common failure mode for
    # stationary recordings with low LF/HF ratio).
    a_p_hat <- exp(0.7 * log(a_p_hat) + 0.3 * log(bc_ap))
    a_s_hat <- exp(0.7 * log(a_s_hat) + 0.3 * log(bc_as))
    if (verbose) message(sprintf(
      "spectral_init: ACF fit (%s) + band centroid ensemble a_p=%.3f a_s=%.3f wp=%.3f",
      biexp$source, a_p_hat, a_s_hat, wp_hat))
  } else {
    # ACF fit failed: fall back entirely to band-centroid estimates
    a_p_hat <- bc_ap;  a_s_hat <- bc_as
    wp_hat  <- bc$free$sigma_p^2 / (2 * bc_ap) /
      (bc$free$sigma_p^2 / (2 * bc_ap) + bc$free$sigma_s^2 / (2 * bc_as))
    ws_hat  <- 1 - wp_hat
    if (verbose) message(sprintf(
      "spectral_init: ACF fit failed; using band centroid a_p=%.3f a_s=%.3f",
      a_p_hat, a_s_hat))
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
