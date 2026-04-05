# ============================================================
# R/model_functions.R
#
# Canonical model infrastructure for the SDE-IG HRV framework.
# Source this file once; it defines all generative-model functions.
# No side-effects: does NOT define PARAMS, INPUT_FN, or RES.
# ============================================================

# ---- Parameter object ----

make_model_params <- function(
    a_p = 2.0, a_s = 0.2,
    sigma_p = 0.30, sigma_s = 0.20,
    mu_0 = 0.85, kappa = 15.0,
    c_p = 0.0, c_s = 0.0) {
  stopifnot(a_p > 0, a_s > 0, a_p > a_s,
            sigma_p > 0, sigma_s > 0,
            mu_0 > 0, kappa > 0)
  list(
    structural = list(a_p = a_p, a_s = a_s,
                      b_ps = 0.0, b_sp = 0.0,
                      k_par = 1.0, k_sym = 1.0),
    free       = list(sigma_p = sigma_p, sigma_s = sigma_s,
                      mu_0 = mu_0, kappa = kappa,
                      c_p = c_p, c_s = c_s)
  )
}

# ---- Log link ----

compute_mu    <- function(p, s, mu_0) mu_0 * exp(p - s)
compute_delta <- function(p, s) s - p

# ---- Exact OU transition kernel ----

ou_exact_step <- function(x_now, a, sigma, c_gain = 0, u_now = 0, dt) {
  e_adt     <- exp(-a * dt)
  mean_next <- x_now * e_adt + (c_gain * u_now / a) * (-expm1(-a * dt))
  var_next  <- (sigma^2 / (2 * a)) * (-expm1(-2 * a * dt))
  mean_next + sqrt(var_next) * rnorm(1)
}

ou_log_transition_density <- function(x_next, x_now, a, sigma,
                                      c_gain = 0, u_now = 0, dt) {
  e_adt     <- exp(-a * dt)
  mean_next <- x_now * e_adt + (c_gain * u_now / a) * (-expm1(-a * dt))
  var_next  <- (sigma^2 / (2 * a)) * (-expm1(-2 * a * dt))
  dnorm(x_next, mean = mean_next, sd = sqrt(var_next), log = TRUE)
}

# ---- Log-space IG computations (Machler 2012) ----

log1mexp <- function(x) {
  stopifnot(all(x < 0))
  ifelse(x >= log(0.5), log(-expm1(x)), log1p(-exp(x)))
}

log_ig_pdf <- function(tau, mu, kappa)
  0.5 * log(kappa) - 0.5 * log(2 * pi * tau^3) -
  kappa * (tau - mu)^2 / (2 * mu^2 * tau)

log_ig_survival <- function(tau, mu, kappa) {
  sq   <- sqrt(kappa / tau)
  z1   <- sq * (tau / mu - 1)
  z2   <- sq * (tau / mu + 1)
  logA <- pnorm(-z1, log.p = TRUE)
  logB <- 2 * kappa / mu + pnorm(-z2, log.p = TRUE)
  d    <- logB - logA
  if (is.nan(d) || d >= 0) return(-Inf)
  logA + log1mexp(d)
}

log_ig_hazard <- function(tau, mu, kappa) {
  if (tau < 1e-9) return(-Inf)
  log_ig_pdf(tau, mu, kappa) - log_ig_survival(tau, mu, kappa)
}

ig_hazard <- function(tau, mu, kappa) exp(log_ig_hazard(tau, mu, kappa))

ig_hazard_peak <- function(mu, kappa) {
  opt <- optimize(function(tau) log_ig_hazard(tau, mu, kappa),
                  interval = c(1e-3, 3 * mu), maximum = TRUE, tol = 1e-6)
  list(tau_peak = opt$maximum, lambda_peak = exp(opt$objective))
}

ig_hazard_bound <- function(mu_min, kappa, safety_factor = 1.1)
  ig_hazard_peak(mu_min, kappa)$lambda_peak * safety_factor

# ---- Survival-ratio event generation ----

p_fire_survival_ratio <- function(tau, mu, kappa, dt) {
  if (tau < 1e-9) return(0)
  ls0 <- log_ig_survival(tau,      mu, kappa)
  ls1 <- log_ig_survival(tau + dt, mu, kappa)
  if (!is.finite(ls0)) return(0)
  min(max(-expm1(ls1 - ls0), 0), 1)
}

lewis_ogata_thinning <- function(time_grid, mu_grid, kappa, lambda_star) {
  T_end   <- tail(time_grid, 1)
  spikes  <- numeric(0)
  t_cur   <- time_grid[1]
  t_last  <- -Inf
  mu_fn   <- approxfun(time_grid, mu_grid, rule = 2)
  while (t_cur < T_end) {
    t_prop  <- t_cur + rexp(1, lambda_star)
    if (t_prop >= T_end) break
    tau_p   <- t_prop - t_last
    lam_p   <- ig_hazard(tau_p, mu_fn(t_prop), kappa)
    if (runif(1) < lam_p / lambda_star) {
      spikes <- c(spikes, t_prop); t_last <- t_prop
    }
    t_cur <- t_prop
  }
  spikes
}

# ---- Observation model interface ----

make_rr_obs_model <- function() {
  list(
    type            = "point_process_ig",
    required_states = c("p", "s"),
    description     = "IG point process; identifiable qty: Delta(t)=s(t)-p(t)",
    log_intensity   = function(tau, p, s, fp) {
      if (tau < 1e-3) return(-Inf)
      log_ig_hazard(tau, compute_mu(p, s, fp$mu_0), fp$kappa)
    }
  )
}

# ---- Main simulation ----
# Returns the simulation list including input_fn for downstream MLE use.

sim_sde_ig <- function(duration, dt, params,
                       input_fn   = function(t) 0,
                       obs_models = list(make_rr_obs_model()),
                       seed       = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- params$structural; fp <- params$free
  n  <- ceiling(duration / dt)
  tg <- seq(0, by = dt, length.out = n)
  p  <- numeric(n); s <- numeric(n)
  mu_v <- numeric(n); dlt <- numeric(n); loglam <- numeric(n)
  spikes <- numeric(0); last_t <- 0
  for (i in seq_len(n - 1)) {
    u <- input_fn(tg[i])
    p[i+1] <- ou_exact_step(p[i], sp$a_p, fp$sigma_p, fp$c_p, u, dt)
    s[i+1] <- ou_exact_step(s[i], sp$a_s, fp$sigma_s, fp$c_s, u, dt)
    mu_v[i+1] <- compute_mu(p[i+1], s[i+1], fp$mu_0)
    dlt[i+1]  <- compute_delta(p[i+1], s[i+1])
    tau       <- tg[i+1] - last_t
    loglam[i+1] <- obs_models[[1]]$log_intensity(tau, p[i+1], s[i+1], fp)
    pf <- p_fire_survival_ratio(tau, mu_v[i+1], fp$kappa, dt)
    if (runif(1) < pf) { spikes <- c(spikes, tg[i+1]); last_t <- tg[i+1] }
  }
  # NOTE: input_fn is stored in the result to enable downstream MLE and
  # identifiability functions to reconstruct the input signal.
  list(time = tg, p = p, s = s, delta = dlt,
       mu = mu_v, lambda = exp(loglam), log_lambda = loglam,
       spikes = spikes, params = params, input_fn = input_fn)
}

# ---- Time-rescaling diagnostics ----

compute_time_rescaling <- function(res) {
  lam <- res$lambda; tg <- res$time; sp <- res$spikes
  dt  <- tg[2] - tg[1]
  lm  <- (lam[-1] + lam[-length(lam)]) / 2
  cml <- c(0, cumsum(lm * dt))
  n   <- length(sp)
  if (n < 2) return(numeric(0))
  vapply(seq_len(n - 1), function(k) {
    i0 <- which.min(abs(tg - sp[k]))
    i1 <- which.min(abs(tg - sp[k + 1]))
    if (i1 <= i0) 0 else cml[i1] - cml[i0]
  }, numeric(1))
}

# ---- Numerical comparison utilities ----

study_firing_bias <- function(mu = 0.85, kappa = 12, dt = 0.005,
                              tau_range = c(0.05, 2.5)) {
  tau_v  <- seq(tau_range[1], tau_range[2], length.out = 600)
  p_surv <- vapply(tau_v, p_fire_survival_ratio, numeric(1),
                   mu = mu, kappa = kappa, dt = dt)
  lam_v  <- vapply(tau_v, ig_hazard, numeric(1), mu = mu, kappa = kappa)
  p_bern <- pmin(lam_v * dt, 1)
  rel_e  <- ifelse(p_surv > 1e-10, (p_bern - p_surv) / p_surv * 100, 0)
  data.frame(tau = tau_v, p_survival = p_surv, p_bernoulli = p_bern,
             lambda = lam_v, rel_error_pct = rel_e)
}

# ============================================================
# Vectorized IG functions
# ============================================================
# The scalar ig_hazard() calls log_ig_survival() which contains
# an early-return branch — not vectorizable.  These V-variants
# implement the same log-space algorithm purely with ifelse / []
# so they work on full simulation grids (n ~ 60 000 – 150 000).
# Used exclusively by the exact conditional MLE in mle.R.

log_ig_pdf_v <- function(tau, mu, kappa) {
  0.5 * log(kappa) - 0.5 * log(2 * pi * tau^3) -
    kappa * (tau - mu)^2 / (2 * mu^2 * tau)
}

log_ig_survival_v <- function(tau, mu, kappa) {
  sq   <- sqrt(kappa / tau)
  z1   <- sq * (tau / mu - 1)
  z2   <- sq * (tau / mu + 1)
  logA <- pnorm(-z1, log.p = TRUE)
  logB <- 2 * kappa / mu + pnorm(-z2, log.p = TRUE)
  d    <- logB - logA
  # log1mexp branch: use the numerically stable variant
  lme  <- ifelse(d >= log(0.5), log(-expm1(d)), log1p(-exp(d)))
  out  <- logA + lme
  out[!is.finite(d) | d >= 0] <- -Inf
  out
}

log_ig_hazard_v <- function(tau, mu, kappa) {
  out        <- log_ig_pdf_v(tau, mu, kappa) -
    log_ig_survival_v(tau, mu, kappa)
  out[tau < 1e-9 | !is.finite(tau)] <- -Inf
  out
}

# ---- Elapsed-time helper ----------------------------------------
# Returns tau(t) = t - t_{last spike before t} for every element of
# time_grid.  Uses findInterval for O(n log m) vectorized computation.

compute_elapsed_times <- function(time_grid, spikes) {
  if (length(spikes) == 0L) return(time_grid - time_grid[1L])
  idx        <- findInterval(time_grid, spikes)  # 0 .. length(spikes)
  last_spike <- c(0, spikes)[idx + 1L]           # 0 before first spike
  time_grid  - last_spike
}
