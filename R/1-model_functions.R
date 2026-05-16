# ============================================================
# R/model_functions.R
#
# Canonical model infrastructure for the SDE-IG HRV framework.
# Source this file once; it defines all generative-model functions.
# No side-effects: does NOT define PARAMS, INPUT_FN, or RES.
# ============================================================

# Minimum physiological a_p / a_s ratio.  Enforced in spectral initialisation,
# Whittle Stage 1, and the full marginal MLE to prevent pole-merging degeneracy.
# Physiological basis: resting HF vagal time-constant (τ_p ≤ 2 s, a_p ≥ 0.5 Hz)
# is at least 3× faster than the sympathetic branch across all known protocols.
SDE_IG_MIN_AP_AS_RATIO <- 3.0

# ---- Double-logistic input function ----
# Smooth, differentiable approximation to a boxcar transient:
#   u(t) ≈ 0  (rest)  →  1  (exercise)  →  0  (recovery)
# with onset at t_on, offset at t_off, and steepness k (Hz).
# At k = 0.1: 10–90% rise ≈ ln(81)/0.1 ≈ 44 s  ≫  τ_p = 0.5 s, τ_s = 5 s.
# From the inference side this looks like naturally slow sympathetic
# accumulation, providing a realistic non-stationarity test without
# an externally visible step. The boxcar is the k → ∞ limit.
#
# Usage:  INPUT_FN <- make_double_logistic(t_on = 300, t_off = 420)

make_double_logistic <- function(t_on, t_off, k = 0.1) {
  force(t_on); force(t_off); force(k)
  function(t) plogis(k * (t - t_on)) * plogis(k * (t_off - t))
}

# ---- Multi-epoch autonomic protocol ----
# Two perturbation epochs with deliberately different time scales:
#   Epoch 1 (t_b_on → t_b_off): SHORT steep burst — primarily drives
#     fast vagal withdrawal, constraining a_p from the rapid onset dynamics.
#   Epoch 2 (t_l_on → t_l_off): SUSTAINED gradual load — primarily drives
#     slow sympathetic accumulation, constraining a_s from the long recovery tail.
#
# This asymmetric design ensures that the marginal log-likelihood surface has
# independent curvature along the a_p and a_s directions: the burst residuals
# carry FIM mass on a_p, while the sustained epoch residuals carry FIM mass on a_s.
# For the conditional MLE a single 60-s window suffices; for the marginal MLE
# the additional a_s epoch substantially reduces the κ–σ_Δ confound.
#
# Default layout (total = 480 s):
#   0 – 120 s : baseline (establishes μ₀, ρ, σ_Δ)
#   120 – 165 s: steep burst (k_b = 0.5 Hz; constrains a_p ≈ 2 Hz)
#   165 – 240 s: first recovery
#   240 – 420 s: sustained load (k_l = 0.1 Hz; constrains a_s ≈ 0.2 Hz)
#   420 – 480 s: second recovery

make_multi_epoch_protocol <- function(
    t_b_on  = 120, t_b_off  = 165, k_b  = 0.5,   # burst
    t_l_on  = 240, t_l_off  = 420, k_l  = 0.1,   # sustained load
    amp_b   = 1.0, amp_l    = 0.7) {              # amplitudes
  force(t_b_on); force(t_b_off); force(k_b)
  force(t_l_on); force(t_l_off); force(k_l)
  force(amp_b);  force(amp_l)
  function(t) {
    burst <- amp_b * plogis(k_b * (t - t_b_on)) * plogis(k_b * (t_b_off - t))
    load  <- amp_l * plogis(k_l * (t - t_l_on))  * plogis(k_l * (t_l_off - t))
    burst + load
  }
}

# ---- Parameter object ----

make_model_params <- function(
    a_p    = 2.0,  a_s   = 0.2,      # self-decay rates (estimated from data)
    a_ps   = 0.0,  a_sp  = 0.0,      # antagonism magnitudes: s→p and p→s inhibition
    sigma_p = 0.30, sigma_s = 0.20,
    mu_0   = 0.85, rho   = 0.20,     # rho = sqrt(mu_0/kappa) = baseline CV of IBI
    c_p    = 0.0,  c_s   = 0.0) {    # input sensitivities (simulation / forced protocols)
  stopifnot(
    a_p > 0, a_s > 0,
    a_ps >= 0, a_sp >= 0,
    a_p * a_s - a_ps * a_sp > 0,     # det(A) > 0: stability
    sigma_p > 0, sigma_s > 0,
    mu_0 > 0, rho > 0
  )
  list(
    structural = list(k_par = 1.0, k_sym = 1.0),  # observation gains; always fixed
    free = list(
      a_p = a_p, a_s = a_s,
      a_ps = a_ps, a_sp = a_sp,
      sigma_p = sigma_p, sigma_s = sigma_s,
      mu_0 = mu_0, rho = rho,
      c_p = c_p, c_s = c_s
    )
  )
}

# Helper: derive kappa from the (mu_0, rho) parameterisation.
# kappa = mu_0 / rho^2  so that  Var(tau_0) = rho^2 * mu_0^2 at baseline.
kappa_from_rho <- function(mu_0, rho) {
  mu_0 / rho^2
}

# ---- Symmetric numerical Hessian via 4-point central differences ----
# Exploits H[i,j] = H[j,i] to halve the number of function evaluations.
# fn   : function of theta returning a scalar (negative log-likelihood)
# theta: parameter vector at which to evaluate
# eps  : step size (default 1e-4 works for log-scale parameterisations)
numerical_hessian <- function(fn, theta, eps = 1e-4) {
  np  <- length(theta)
  H   <- matrix(0, np, np)
  ep  <- diag(eps, np)              # column i is the i-th canonical basis vector × eps
  for (i in seq_len(np)) {
    for (j in i:np) {
      H[i, j] <- (fn(theta + ep[, i] + ep[, j]) - fn(theta + ep[, i] - ep[, j]) -
                    fn(theta - ep[, i] + ep[, j]) + fn(theta - ep[, i] - ep[, j])) /
                   (4 * eps^2)
      H[j, i] <- H[i, j]
    }
  }
  H
}

# ---- Log link ----

compute_mu    <- function(p, s, mu_0) {
  mu_0 * exp(p - s)
}

compute_delta <- function(p, s) {
  s - p
}

# ---- 2×2 analytical matrix exponential (no external dependencies) ----
# Uses the Cayley-Hamilton spectral decomposition:
#   e^{At} = α₀(t)I + α₁(t)(A − (tr/2)I)
# covering three eigenvalue regimes: real distinct, repeated, complex.

mat2x2_exp <- function(A, t) {
  tr_A <- A[1L, 1L] + A[2L, 2L]
  det_A <- A[1L, 1L] * A[2L, 2L] - A[1L, 2L] * A[2L, 1L]
  disc  <- tr_A^2 - 4 * det_A      # (λ1 − λ2)²

  ht  <- tr_A * t / 2
  eht <- exp(ht)
  B   <- A - (tr_A / 2) * diag(2L) # A − (tr/2)I

  tol <- 1e-10 * (abs(tr_A)^2 + 1)
  if (abs(disc) <= tol) {
    eht * (diag(2L) + B * t)
  } else if (disc > 0) {
    sq   <- sqrt(disc) / 2
    # Eigenvalue form: avoids cosh/sinh overflow × eht underflow → NaN.
    # e^{At} = e^{λ₁t}/2 (I + B/sq) + e^{λ₂t}/2 (I - B/sq)
    # λ₁ = tr_A/2 + sq < 0,  λ₂ = tr_A/2 − sq < 0  (both negative by stability).
    e1  <- exp((ht + sq * t) - log(2))   # exp(λ₁ t) / 2
    e2  <- exp((ht - sq * t) - log(2))   # exp(λ₂ t) / 2
    Bsq <- B / sq
    e1 * (diag(2L) + Bsq) + e2 * (diag(2L) - Bsq)
  } else {
    om <- sqrt(-disc) / 2
    eht * (cos(om * t) * diag(2L) + sin(om * t) / om * B)
  }
}

# ---- Exact stationary covariance for the 2×2 coupled OU system ----
# Analytic solution of A P∞ + P∞ Aᵀ + ΣΣᵀ = 0 (algebraic Lyapunov equation).
# Reduces to diag(σ_p²/2a_p, σ_s²/2a_s) when a_ps = a_sp = 0.
ou_stationary_cov <- function(a_p, a_s, a_ps, a_sp, sigma_p, sigma_s) {
  D  <- a_p * a_s - a_ps * a_sp    # det(−A) > 0 by stability
  Tr <- a_p + a_s                  # tr(−A) > 0
  p11 <- ((a_s^2 + D) * sigma_p^2 + a_ps^2 * sigma_s^2) / (2 * D * Tr)
  p22 <- ((a_p^2 + D) * sigma_s^2 + a_sp^2 * sigma_p^2) / (2 * D * Tr)
  p12 <- -(a_ps * p22 + a_sp * p11) / Tr
  matrix(c(p11, p12, p12, p22), 2L, 2L)
}

# ---- Pre-compute all objects needed for simulation and likelihood ----
# Calling this once per parameter vector avoids redundant matrix exponentials.

ou_coupled_matrices <- function(a_p, a_s, a_ps, a_sp,
                                sigma_p, sigma_s, dt) {
  # A = [[-a_p, -a_ps], [-a_sp, -a_s]]:
  #   dp = -a_p*p - a_ps*s  (s inhibits p — accentuated antagonism)
  #   ds = -a_sp*p - a_s*s  (p inhibits s — accentuated antagonism)
  # Stability: det(A) = a_p*a_s - a_ps*a_sp > 0 (enforced by make_model_params)
  A_mat <- matrix(c(-a_p, -a_sp, -a_ps, -a_s), 2L, 2L)
  F_mat <- mat2x2_exp(A_mat, dt)
  # Reuse F_mat: Q(dt) = P_inf - F(dt) P_inf F(dt)^T (exact Lyapunov identity).
  # Avoids the redundant mat2x2_exp call that ou_coupled_Q would otherwise make.
  P_inf <- ou_stationary_cov(a_p, a_s, a_ps, a_sp, sigma_p, sigma_s)
  Q_raw <- P_inf - F_mat %*% P_inf %*% t(F_mat)
  Q_mat <- 0.5 * (Q_raw + t(Q_raw))

  # 2×2 symmetric PSD check via Sylvester's criterion: trace > 0 and det > 0.
  # Avoids eigen() entirely — O(1) arithmetic replaces iterative QR.
  det_Q <- Q_mat[1L, 1L] * Q_mat[2L, 2L] - Q_mat[1L, 2L]^2
  if (Q_mat[1L, 1L] <= 0 || det_Q <= 0) {
    tr_half <- (Q_mat[1L, 1L] + Q_mat[2L, 2L]) / 2
    min_eig <- tr_half - sqrt(max(tr_half^2 - det_Q, 0))
    shift   <- abs(min_eig) + 1e-10
    Q_mat   <- Q_mat + shift * diag(2L)
    det_Q   <- Q_mat[1L, 1L] * Q_mat[2L, 2L] - Q_mat[1L, 2L]^2
  }

  # Closed-form 2×2 inverse and log-determinant: exact, no iteration.
  log_det_Q <- log(det_Q)

  Q_inv     <- matrix(c( Q_mat[2L, 2L], -Q_mat[1L, 2L],
                         -Q_mat[1L, 2L],  Q_mat[1L, 1L]), 2L, 2L) / det_Q

  Q_chol    <- tryCatch(chol(Q_mat),
                        error = function(e) chol(Q_mat + 1e-9 * diag(2L)))

  D_mat <- a_p * a_s - a_ps * a_sp

  A_inv <- matrix(c(-a_s, a_sp, a_ps, -a_p), 2L, 2L) / D_mat

  list(A = A_mat, F = F_mat, Ft = t(F_mat), FmI = F_mat - diag(2L),
       Q = Q_mat, Q_inv = Q_inv, log_det_Q = log_det_Q,
       Q_chol = Q_chol, A_inv = A_inv)
}

# ---- Joint 2D OU simulation step ----

ou_coupled_step <- function(x_now, mats, c_vec, u_now) {
  mean_next <- as.vector(mats$F %*% x_now)
  if (u_now != 0) {
    d         <- as.vector(mats$FmI %*% mats$A_inv %*% (c_vec * u_now))
    mean_next <- mean_next + d
  }
  mean_next + as.vector(t(mats$Q_chol) %*% rnorm(2L))
}

# ---- Exact OU transition kernel ----

ou_exact_step <- function(x_now, a, sigma, c_gain = 0, u_now = 0, dt) {
  e_adt     <- exp(-a * dt)
  mean_next <- x_now * e_adt + (c_gain * u_now / a) * (-expm1(-a * dt))
  var_next  <- (sigma^2 / (2 * a)) * (-expm1(-2 * a * dt))
  mean_next + sqrt(var_next) * rnorm(1)
}

# ---- Log-space IG computations (Machler 2012) ----

log1mexp <- function(x) {
  # log(1 - exp(x)), numerically stable for x < 0.
  # Values >= 0 indicate upstream overflow; clamp element-wise and warn.
  bad <- is.finite(x) & x >= 0
  if (any(bad, na.rm = TRUE)) {
    warning("log1mexp: ", sum(bad), " value(s) >= 0; clamping to -eps.")
    x <- ifelse(bad, -.Machine$double.eps, x)
  }
  ifelse(x >= log(0.5), log(-expm1(x)), log1p(-exp(x)))
}

log_ig_pdf <- function(tau, mu, kappa) {
  ifelse(
    is.finite(tau) & tau > 0,
    0.5 * log(kappa) - 0.5 * log(2 * pi * tau^3) -
      kappa * (tau - mu)^2 / (2 * mu^2 * tau),
    -Inf
  )
}

log_ig_survival <- function(tau, mu, kappa) {
  # Vectorized over tau and mu (kappa remains scalar).
  n   <- max(length(tau), length(mu))
  tau <- rep_len(tau, n)
  mu  <- rep_len(mu,  n)
  out <- rep.int(0.0, n)      # log S = 0 (S = 1) for tau <= 0 or non-finite mu

  pos <- is.finite(tau) & tau > 0 & is.finite(mu)
  if (!any(pos)) return(out)

  tp  <- tau[pos];  mp <- mu[pos]
  sq  <- sqrt(kappa / tp)
  z1  <- sq * (tp / mp - 1)
  z2  <- sq * (tp / mp + 1)

  s   <- numeric(sum(pos))
  bel <- tp < mp               # --- branch: tau < mu, compute F then S ---

  if (any(bel)) {
    lF1 <- pnorm( z1[bel], log.p = TRUE)
    lF2 <- 2 * kappa / mp[bel] + pnorm(-z2[bel], log.p = TRUE)
    lF  <- ifelse(lF1 >= lF2,
                  lF1 + log1p(exp(lF2 - lF1)),
                  lF2 + log1p(exp(lF1 - lF2)))
    bad <- !is.finite(lF) | lF >= 0
    if (any(bad)) {
      warning(sprintf(
        "log_ig_survival: %d value(s) with log_F >= 0 (first: tau=%.3g, mu=%.3g); clamping.",
        sum(bad), tp[bel][bad][1L], mp[bel][bad][1L]))
      lF[bad] <- -.Machine$double.eps
    }
    s[bel] <- log1mexp(lF)
  }

  abv <- !bel                  # --- branch: tau >= mu, compute S directly ---
  if (any(abv)) {
    lA     <- pnorm(-z1[abv], log.p = TRUE)
    d      <- 2 * kappa / mp[abv] + pnorm(-z2[abv], log.p = TRUE) - lA
    s[abv] <- ifelse(is.nan(d) | d >= 0, -Inf, lA + log1mexp(d))
  }

  out[pos] <- s
  out
}

log_ig_hazard <- function(tau, mu, kappa) {
  lh <- log_ig_pdf(tau, mu, kappa) - log_ig_survival(tau, mu, kappa)
  lh[!is.finite(tau) | tau < 1e-9] <- -Inf
  lh
}

ig_hazard <- function(tau, mu, kappa) exp(log_ig_hazard(tau, mu, kappa))

ig_hazard_peak <- function(mu, kappa) {
  opt <- optimize(function(tau) log_ig_hazard(tau, mu, kappa),
                  interval = c(1e-3, 3 * mu), maximum = TRUE, tol = 1e-6)
  list(tau_peak = opt$maximum, lambda_peak = exp(opt$objective))
}

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

# ---- Observation model interface ----

make_rr_obs_model <- function() {
  list(
    type            = "point_process_ig",
    required_states = c("p", "s"),
    description     = "IG point process; identifiable qty: Delta(t)=s(t)-p(t)",
    log_intensity   = function(tau, p, s, fp) {
      if (tau < 1e-3) return(-Inf)
      log_ig_hazard(tau, compute_mu(p, s, fp$mu_0), kappa_from_rho(fp$mu_0, fp$rho))
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

  # Initialize at the stationary mean of each OU process.
  # When input is non-zero at t=0 the equilibrium differs, but the
  # stationary variance is sigma^2/(2a) regardless of input.
  # Drawing from the marginal N(0, sigma^2/(2a)) avoids a systematic
  # transient that would bias short simulations.
  p  <- numeric(n)
  s  <- numeric(n)


  # Draw from the exact bivariate stationary distribution.
  # Uses the closed-form P∞; reduces to independent draws when a_ps = a_sp = 0.
  P_stat <- ou_stationary_cov(fp$a_p, fp$a_s, fp$a_ps, fp$a_sp,
                              fp$sigma_p, fp$sigma_s)
  x0 <- tryCatch(
    as.vector(t(chol(P_stat)) %*% rnorm(2L)),
    error = function(e) c(rnorm(1L, 0, sqrt(fp$sigma_p^2 / (2 * fp$a_p))),
                          rnorm(1L, 0, sqrt(fp$sigma_s^2 / (2 * fp$a_s))))
  )
  p[1L] <- x0[1L];  s[1L] <- x0[2L]

  # Pre-compute coupled transition matrices once if coupling is non-zero.
  # Falls back to the scalar exact-step fast path when b_ps = b_sp = 0.
  use_coupled <- (abs(fp$a_ps) + abs(fp$a_sp)) > 1e-12
  if (use_coupled) {
    cmats <- ou_coupled_matrices(fp$a_p, fp$a_s, fp$a_ps, fp$a_sp,
                                 fp$sigma_p, fp$sigma_s, dt)
    c_vec_ou <- c(fp$c_p, fp$c_s)
  }

  mu_v <- numeric(n); dlt <- numeric(n); loglam <- numeric(n)
  loglam[1L] <- -Inf
  mu_v[1L]   <- compute_mu(p[1L], s[1L], fp$mu_0)
  dlt[1L]    <- compute_delta(p[1L], s[1L])

  max_spikes <- n
  spikes     <- numeric(max_spikes)
  n_spikes   <- 0L
  last_t     <- 0

  kap      <- kappa_from_rho(fp$mu_0, fp$rho)
  # Pre-evaluate input on the full time grid; input_fn (make_double_logistic,
  # make_multi_epoch_protocol) is vectorized via plogis. Avoids n-1 closure calls.
  u_sim <- input_fn(tg[-length(tg)])
  stopifnot(length(u_sim) == length(tg) - 1L)

  for (i in seq_len(n - 1)) {
    u <- u_sim[i]
    tau_start <- tg[i] - last_t
    tau_end   <- tg[i+1] - last_t

    # Firing probability evaluated at START of interval [tg[i], tg[i+1]]
    pf <- p_fire_survival_ratio(tau_start, tau_end, mu_v[i], kap)
    if (runif(1) < pf) {
      n_spikes <- n_spikes + 1L
      spikes[n_spikes] <- tg[i+1]
      last_t <- tg[i+1]
    }

    # State advance AFTER firing decision
    if (use_coupled) {
      xnew   <- ou_coupled_step(c(p[i], s[i]), cmats, c_vec_ou, u)
      p[i+1] <- xnew[1L];  s[i+1] <- xnew[2L]
    } else {
      p[i+1] <- ou_exact_step(p[i], fp$a_p, fp$sigma_p, fp$c_p, u, dt)
      s[i+1] <- ou_exact_step(s[i], fp$a_s, fp$sigma_s, fp$c_s, u, dt)
    }
    mu_v[i+1] <- compute_mu(p[i+1], s[i+1], fp$mu_0)
    dlt[i+1]  <- compute_delta(p[i+1], s[i+1])
    loglam[i+1] <- obs_models[[1]]$log_intensity(tg[i+1] - last_t, p[i+1], s[i+1], fp)
  }
  spikes <- spikes[seq_len(n_spikes)]
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
  fp    <- res$params$free
  kappa <- kappa_from_rho(fp$mu_0, fp$rho)
  n     <- length(sp)
  if (n < 2L) return(numeric(0L))

  i_bounds <- findInterval(sp, tg)   # length n, O(n log m) total

  # Lengths of each IBI's grid-point span
  lens <- i_bounds[-1L] - i_bounds[-n]        # i1 - i0 for each IBI
  bad  <- lens <= 0L
  if (any(bad)) {
    warning(sprintf(
      "compute_time_rescaling: %d IBI(s) shorter than dt; returning NA for those.",
      sum(bad)))
  }

  # Build flat index vectors for all valid IBIs at once
  valid     <- which(!bad)
  all_i1    <- i_bounds[-1L][valid]
  all_i0    <- i_bounds[-n][valid]
  all_lens  <- lens[valid]

  # sequence() builds (i0+1):i1 for every IBI without an R loop
  flat_end  <- sequence(all_lens, from = all_i0 + 1L)
  flat_start <- flat_end - 1L
  ibi_id    <- rep(valid, times = all_lens)

  # Vectorised spike-relative times
  tau_ends   <- tg[flat_end]   - sp[ibi_id]
  tau_starts <- tg[flat_start] - sp[ibi_id]
  # Clamp the very first entry of each IBI (tg[i0] <= sp[k] so tau_start[1] <= 0)
  first_of_ibi              <- !duplicated(ibi_id)
  tau_starts[first_of_ibi]  <- 0

  mu_j <- mu_v[flat_start]

  delta_ll <- log_ig_survival(tau_starts, mu_j, kappa) -
    log_ig_survival(tau_ends,   mu_j, kappa)

  # Aggregate per-IBI
  Lk        <- rep(NA_real_, n - 1L)
  Lk[valid] <- as.numeric(rowsum(delta_ll, ibi_id)[, 1L])
  Lk
}

# ---- Numerical comparison utilities ----

study_firing_bias <- function(mu = 0.85, kappa = 12, dt = 0.005,
                              tau_range = c(0.05, 2.5)) {
  tau_v  <- seq(tau_range[1], tau_range[2], length.out = 600)
  p_surv <- vapply(tau_v, function(t) p_fire_survival_ratio(t, t + dt, mu, kappa), numeric(1))

  # Evaluate hazard at step entry (tau_v) — consistent with the survival-ratio simulation.
  # Compares p_bern = lambda(tau_start)*dt  against  p_surv = 1 - S(tau+dt)/S(tau).
  lam_start <- vapply(tau_v, ig_hazard, numeric(1), mu = mu, kappa = kappa)
  p_bern    <- pmin(lam_start * dt, 1)

  rel_e    <- ifelse(p_surv > 1e-10, (p_bern - p_surv) / p_surv * 100, 0)
  lam_plot <- lam_start   # reuse; already evaluated at tau_v

  data.frame(tau = tau_v, p_survival = p_surv, p_bernoulli = p_bern,
             lambda = lam_plot, rel_error_pct = rel_e)
}
