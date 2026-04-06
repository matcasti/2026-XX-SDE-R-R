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
  if (tau <= 0) return(0)
  sq   <- sqrt(kappa / tau)
  z1   <- sq * (tau / mu - 1)
  z2   <- sq * (tau / mu + 1)

  if (tau < mu) {
    logF1 <- pnorm( z1, log.p = TRUE)           # log Φ(z1),  z1 < 0
    logF2 <- 2 * kappa / mu + pnorm(-z2, log.p = TRUE)  # log[exp(2κ/μ)Φ(−z2)]
    # log F(τ) = log[exp(logF1) + exp(logF2)] via log-sum-exp
    log_F  <- if (logF1 >= logF2)
                logF1 + log1p(exp(logF2 - logF1))
              else
                logF2 + log1p(exp(logF1 - logF2))
    # log S(τ) = log(1 − F(τ))

    if (log_F >= 0) return(-Inf)   # numerical edge case: F ≥ 1

    return(log1mexp(log_F))    # log1mexp(x) = log(1−eˣ), x < 0 required
  } else {
    logA <- pnorm(-z1, log.p = TRUE)
    logB <- 2 * kappa / mu + pnorm(-z2, log.p = TRUE)
    d    <- logB - logA
    if (is.nan(d) || d >= 0) return(-Inf)
    return(logA + log1mexp(d))
  }
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

p_fire_survival_ratio <- function(tau_start, tau_end, mu, kappa) {
  ls0 <- log_ig_survival(tau_start, mu, kappa)
  # S(tau_start) ≈ 0: event is overdue, fire with certainty
  if (!is.finite(ls0)) return(1)
  ls1 <- log_ig_survival(tau_end, mu, kappa)
  # S(tau_end) ≈ 0 but S(tau_start) > 0: certain fire in this step
  if (!is.finite(ls1)) return(1)
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

  loglam[1L] <- -Inf   # hazard at tau=0 is zero by IG first-passage construction

  spikes <- numeric(0); last_t <- 0
  for (i in seq_len(n - 1)) {
    u <- input_fn(tg[i])
    p[i+1] <- ou_exact_step(p[i], sp$a_p, fp$sigma_p, fp$c_p, u, dt)
    s[i+1] <- ou_exact_step(s[i], sp$a_s, fp$sigma_s, fp$c_s, u, dt)
    mu_v[i+1] <- compute_mu(p[i+1], s[i+1], fp$mu_0)
    dlt[i+1]  <- compute_delta(p[i+1], s[i+1])

    tau_start <- tg[i] - last_t
    tau_end   <- tg[i+1] - last_t

    loglam[i+1] <- obs_models[[1]]$log_intensity(tau_end, p[i+1], s[i+1], fp)

    pf <- p_fire_survival_ratio(tau_start, tau_end, mu_v[i+1], fp$kappa)
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
  sp    <- res$spikes
  tg    <- res$time
  dt    <- tg[2L] - tg[1L]
  mu_v  <- res$mu
  kappa <- res$params$free$kappa
  n     <- length(sp)
  if (n < 2L) return(numeric(0L))

  vapply(seq_len(n - 1L), function(k) {
    i0 <- which.min(abs(tg - sp[k]))
    i1 <- which.min(abs(tg - sp[k + 1L]))
    if (i1 <= i0) return(0)

    # Grid points strictly inside (sp[k], sp[k+1]], i.e. the END of each dt step
    idx      <- (i0 + 1L):i1
    tau_ends <- tg[idx] - sp[k]    # elapsed time at end of each step

    # ΔΛ_j = log S(τ_start; μ_j, κ) − log S(τ_end; μ_j, κ)  ≥ 0
    # Mirrors exactly: P_j = 1 − S(τ_end)/S(τ_start) used in sim_sde_ig().
    # At j=1: τ_start = 0, log_ig_survival(0, ...) = 0 (S(0) = 1) ✓
    sum(vapply(seq_along(idx), function(j) {
      tau_end   <- tau_ends[j]
      tau_start <- max(tau_end - dt, 0)
      mu_j      <- mu_v[idx[j]]
      log_ig_survival(tau_start, mu_j, kappa) -
        log_ig_survival(tau_end,   mu_j, kappa)
    }, numeric(1L)))
  }, numeric(1L))
}

# ---- Numerical comparison utilities ----

study_firing_bias <- function(mu = 0.85, kappa = 12, dt = 0.005,
                              tau_range = c(0.05, 2.5)) {
  tau_v  <- seq(tau_range[1], tau_range[2], length.out = 600)
  p_surv <- vapply(tau_v, function(t) p_fire_survival_ratio(t, t + dt, mu, kappa), numeric(1))

  # Evaluate hazard at midpoint to isolate the Bernoulli truncation bias
  lam_start <- vapply(tau_v, ig_hazard, numeric(1), mu = mu, kappa = kappa)
  p_bern    <- pmin(lam_start * dt, 1)

  rel_e    <- ifelse(p_surv > 1e-10, (p_bern - p_surv) / p_surv * 100, 0)
  lam_plot <- lam_start   # reuse; already evaluated at tau_v

  data.frame(tau = tau_v, p_survival = p_surv, p_bernoulli = p_bern,
             lambda = lam_plot, rel_error_pct = rel_e)
}
