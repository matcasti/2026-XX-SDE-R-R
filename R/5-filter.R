# ============================================================
# R/filter.R
#
# Spike-to-spike Unscented Kalman Filter for the SDE-IG model.
# Requires: model_functions.R and spectral_init.R sourced first.
#
# ALGORITHM
# ---------
# For each IBI τ_k = t_k − t_{k−1}:
#
#  1. PREDICT (t_{k-1} → t_k)
#     Exact OU Gaussian transition — linear, so no sigma-point
#     approximation needed at this step.
#     m_{k|k-1}, P_{k|k-1} in closed form.
#
#  2. UPDATE (observation y_k = τ_k at spike t_k)
#     h(x) = μ₀ exp(p − s) is nonlinear → UKF sigma-point quadrature.
#     Measurement noise R = μ̂³/κ (IG variance outer-linearized at the
#     predicted mean, so R adapts with the autonomic state).
#     Marginal LL contribution: log N(τ_k; μ̂_k, S_k).
#
# IDENTIFIABILITY GEOMETRY
# -------------------------
# Only Δ(t) = s − p is identified by a single RR stream.  After each
# update the posterior variance along e₋ = [−1, 1]ᵀ/√2 (the Δ direction)
# contracts; variance along e₊ = [1, 1]ᵀ/√2 (common mode) remains large,
# governed by the OU prior.  This geometry emerges automatically from the
# filter — no constraint is imposed by hand.
# ============================================================

# ---- UKF tuning (Wan & van der Merwe 2000) ----
# L = 2 (state dim: p, s).  alpha = 0.1 gives moderate sigma-point spread;
# beta = 2 is optimal for Gaussian priors; kappa = 0 is standard.

.ukf_weights <- local({
  # Symmetric cubature rule: alpha=1, kappa=0, beta=2 (Wan & van der Merwe 2000).
  # Gives c = L = 2, no negative weights (Wm[0] = 0, Wm[i] = 0.25),
  # and Wc[0] = 2 (beta correction for Gaussian kurtosis).
  # Sigma-point spread = sqrt(L * P), providing better nonlinear coverage than
  # the alpha=0.1 rule (spread = sqrt(0.02 * P)) while avoiding Wm[0] = -99.
  L <- 2L; alpha <- 1.0; beta <- 2.0; kappa <- 0.0
  lambda  <- alpha^2 * (L + kappa) - L      # = 1*2 - 2 = 0
  c_val   <- L + lambda                     # = 2
  n_pts   <- 2L * L + 1L
  Wm <- c(lambda / c_val, rep(0.5 / c_val, 2L * L))   # [0, 0.25, 0.25, 0.25, 0.25]
  Wc <- c(lambda / c_val + (1 - alpha^2 + beta), rep(0.5 / c_val, 2L * L))  # [2, 0.25, ...]
  list(lambda = lambda, c = c_val, Wm = Wm, Wc = Wc, L = L, n_pts = n_pts)
})

# ---- Sigma-point generation (scaled, Cholesky-based) ----

# L = 2 is structural (fixed by .ukf_weights); unroll the loop entirely.
.sigma_points <- function(m, P) {
  A <- tryCatch(
    t(chol(.ukf_weights$c * P)),
    error = function(e) t(chol(.ukf_weights$c * (P + 1e-8 * diag(2L))))
  )
  cbind(m, m + A[, 1L], m + A[, 2L], m - A[, 1L], m - A[, 2L])
}

# ---- Exact OU prediction ----
# The OU transition is linear-Gaussian, so the prediction step is an
# exact Kalman prediction — no approximation.
#
# Returns: list(m, P)

.ou_predict <- function(m0, P0, tau, fp, u = 0) {
  a_p <- fp$a_p;  a_s <- fp$a_s

  if ((abs(fp$a_ps) + abs(fp$a_sp)) > 1e-12) {
    A_mat <- matrix(c(-a_p, -fp$a_sp, -fp$a_ps, -a_s), 2L, 2L)
    F_mat <- mat2x2_exp(A_mat, tau)
    # Reuse F_mat already computed; derive Q via the Lyapunov identity
    # to avoid the redundant mat2x2_exp call inside ou_coupled_Q.
    P_inf <- ou_stationary_cov(a_p, fp$a_s, fp$a_ps, fp$a_sp,
                               fp$sigma_p, fp$sigma_s)
    Q_raw <- P_inf - F_mat %*% P_inf %*% t(F_mat)
    Q_mat <- 0.5 * (Q_raw + t(Q_raw))

    det_Q <- Q_mat[1L, 1L] * Q_mat[2L, 2L] - Q_mat[1L, 2L]^2
    if (Q_mat[1L, 1L] <= 0 || det_Q <= 0) {
      tr_h  <- (Q_mat[1L, 1L] + Q_mat[2L, 2L]) / 2
      shift <- abs(tr_h - sqrt(max(tr_h^2 - det_Q, 0))) + 1e-10
      Q_mat <- Q_mat + shift * diag(2L)
    }

    d_vec <- c(0, 0)
    if (u != 0) {
      D_loc <- fp$a_p * fp$a_s - fp$a_ps * fp$a_sp
      A_inv <- matrix(c(-a_s, fp$a_sp, fp$a_ps, -a_p), 2L, 2L) / D_loc
      d_vec <- as.vector((F_mat - diag(2L)) %*% A_inv %*% c(fp$c_p * u, fp$c_s * u))
    }

    m1 <- as.vector(F_mat %*% m0) + d_vec
    P1 <- F_mat %*% P0 %*% t(F_mat) + Q_mat
    P1 <- 0.5 * (P1 + t(P1))
    return(list(m = m1, P = P1))
  }

  # Uncoupled fast path: exact scalar formulas (no matrix operations).
  F_p <- exp(-a_p * tau);  F_s <- exp(-a_s * tau)
  d_p <- (fp$c_p * u / a_p) * (-expm1(-a_p * tau))
  d_s <- (fp$c_s * u / a_s) * (-expm1(-a_s * tau))
  m1  <- c(m0[1L] * F_p + d_p, m0[2L] * F_s + d_s)

  Q_p <- (fp$sigma_p^2 / (2 * a_p)) * (-expm1(-2 * a_p * tau))
  Q_s <- (fp$sigma_s^2 / (2 * a_s)) * (-expm1(-2 * a_s * tau))
  P1  <- matrix(0, 2L, 2L)
  P1[1L, 1L] <- F_p^2 * P0[1L, 1L] + Q_p
  P1[2L, 2L] <- F_s^2 * P0[2L, 2L] + Q_s
  P1[1L, 2L] <- F_p * F_s * P0[1L, 2L]
  P1[2L, 1L] <- P1[1L, 2L]
  list(m = m1, P = P1)
}

# ---- UKF update at a spike ----
# Observation: τ_k ~ IG(μ₀ exp(p − s), κ).
# Observation function: h(x) = μ₀ exp(p − s).
# Noise variance: R = μ̂³/κ (IG variance outer-linearized at predicted mean).
#
# Returns: list(m, P, innov, S, mu_hat, ll)

.ukf_update <- function(m1, P1, tau_k, fp, kappa) {
  mu_0  <- fp$mu_0

  pts    <- .sigma_points(m1, P1)
  h_pts  <- mu_0 * exp(pts[1L, ] - pts[2L, ])   # h(sigma_i)

  mu_hat <- sum(.ukf_weights$Wm * h_pts)
  if (!is.finite(mu_hat) || mu_hat <= 0) {
    mu_det <- mu_0 * exp(m1[1L] - m1[2L])
    return(list(m = m1, P = P1, innov = 0, S = 1,
                mu_hat = if (is.finite(mu_det) && mu_det > 0) mu_det else mu_0,
                ll = -Inf))
  }
  R <- max(mu_hat^3 / kappa, 1e-6)

  dh   <- h_pts - mu_hat
  S    <- sum(.ukf_weights$Wc * dh^2) + R

  # Cross-covariance state–observation (2×1)
  dx   <- pts - m1                                 # 2 × n_pts
  C    <- as.vector(dx %*% (.ukf_weights$Wc * dh))  # 2×n_pts · n_pts = 2×1

  K      <- C / S
  innov  <- tau_k - mu_hat

  m2     <- m1 + K * innov
  P2_raw <- P1 - outer(K, K) * S
  P2     <- 0.5 * (P2_raw + t(P2_raw)) + 1e-9 * diag(.ukf_weights$L)

  # Wm[1] = 0 (central point excluded); Wm[2:5] = 0.25 (equal) → log-mean of 4 values.
  ll_sigma <- log_ig_pdf(tau_k, h_pts[-1L], kappa)   # 4 off-centre sigma-point LLs
  finite   <- is.finite(ll_sigma)
  ll_k     <- if (any(finite)) {
    lmax <- max(ll_sigma[finite])
    lmax + log(sum(exp(ll_sigma[finite] - lmax)) / 4L)
  } else {
    -Inf
  }

  list(m = m2, P = P2, innov = innov, S = S, mu_hat = mu_hat, ll = ll_k)
}

# ---- Full spike-to-spike filter ----
#
# Arguments
#   spikes   : ascending numeric vector of spike times
#   params   : make_model_params() output
#   input_fn : function(t) → scalar; exogenous drive at time t
#   m0, P0   : optional initial state mean (length-2) and covariance (2×2)
#              Defaults: stationary OU prior (zero mean, variance σ²/(2a))
#
# Returns
#   m_filt     : (N−1)×2 matrix of filtered state means [p, s]
#   P_filt     : list of (N−1) filtered covariance matrices
#   delta_filt : filtered Δ̂(t) = ŝ − p̂ at each beat time
#   mu_filt    : filtered mean IBI μ̂_k at each beat time
#   beat_times : spike times t_2 … t_N (the IBI observation times)
#   innov      : innovation sequence τ_k − μ̂_k
#   S_innov    : innovation variances
#   ll         : total marginal log-likelihood (sum over beats)
#   ll_vec     : per-beat marginal log-likelihood contributions

pp_ukf <- function(spikes,
                   params,
                   input_fn = function(t) 0,
                   m0       = NULL,
                   P0       = NULL) {
  sp <- params$structural
  fp <- params$free
  kap_const <- kappa_from_rho(fp$mu_0, fp$rho)
  N  <- length(spikes)
  if (N < 2L) stop("pp_ukf: need at least 2 spikes.")

  if (is.null(m0)) m0 <- c(0, 0)
  if (is.null(P0)) P0 <- ou_stationary_cov(fp$a_p, fp$a_s, fp$a_ps, fp$a_sp,
                                           fp$sigma_p, fp$sigma_s)

  n_ibi      <- N - 1L
  m_filt     <- matrix(0, n_ibi, 2L, dimnames = list(NULL, c("p", "s")))
  P_filt     <- vector("list", n_ibi)
  delta_filt <- mu_filt <- innov_v <- S_v <- ll_v <- numeric(n_ibi)

  m_cur <- m0; P_cur <- P0

  for (k in seq_len(n_ibi)) {
    tau_k <- spikes[k + 1L] - spikes[k]
    # Evaluate input at the START of the IBI; treat as constant over [t_{k-1}, t_k].
    # For the double-logistic protocol (10--90% rise ≈ 44 s) and IBIs ~0.85 s,
    # the resulting approximation error in the input-driven mean shift d(τ_k)
    # is O(|u̇|·τ_k²/2), which is negligible everywhere except at the onset/offset
    # transition (≈ 0.5% error at the steepest point for the reference parameters).
    u_k   <- input_fn(spikes[k])

    pred <- .ou_predict(m_cur, P_cur, tau_k, fp, u = u_k)
    upd  <- .ukf_update(pred$m, pred$P, tau_k, fp, kap_const)

    m_filt[k, ]   <- upd$m
    P_filt[[k]]   <- upd$P
    delta_filt[k] <- upd$m[2L] - upd$m[1L]   # ŝ − p̂
    mu_filt[k]    <- upd$mu_hat
    innov_v[k]    <- upd$innov
    S_v[k]        <- upd$S
    ll_v[k]       <- upd$ll

    m_cur <- upd$m; P_cur <- upd$P
  }

  list(
    m_filt     = m_filt,
    P_filt     = P_filt,
    delta_filt = delta_filt,
    mu_filt    = mu_filt,
    beat_times = spikes[-1L],
    innov      = innov_v,
    S_innov    = S_v,
    ll         = sum(ll_v),
    ll_vec     = ll_v,
    n_beats    = n_ibi,
    params     = params
  )
}

# ---- Marginal MLE over static parameters ----
#
# Maximizes the filter log-likelihood over the free parameter set.
# Optimization is unconstrained in log-space for positive parameters.
#
# Arguments
#   spikes         : spike time vector
#   params_init    : starting-value parameter object (make_model_params())
#   input_fn       : exogenous drive function
#   optimize_gains : if TRUE, also optimize c_p and c_s (requires non-trivial input)
#   verbose        : print L-BFGS-B progress
#
# Returns
#   params_hat  : estimated parameter object
#   free_hat    : named list of estimated free parameters
#   ll          : marginal log-likelihood at MLE
#   convergence : optim() convergence code (0 = success)
#   filter      : final filter run at the MLE

pp_mle <- function(spikes,
                   params_init,
                   input_fn = function(t) 0,
                   verbose  = FALSE) {

  # ---- Data-driven initialisation ----
  # spectral_init() estimates starting values from the tachogram's
  # second-order structure (biexponential ACF, SDNN, RMSSD, mean).
  # This replaces hard-coded literature defaults and substantially
  # reduces dependence of L-BFGS-B on the choice of params_init.
  # Coupling terms (a_ps, a_sp) and input gains (c_p, c_s) are
  # carried over from params_init unchanged, since spectral_init
  # does not estimate them.
  rr_for_init <- diff(spikes)

  # Summary statistics fixed from data — captured by the pack/unpack closures.
  # mu_bar: sample mean IBI (seconds).  cv_obs: observed coefficient of variation
  # of IBIs (floor at 0.01 prevents division by zero on pathologically regular data).

  mu_bar_raw <- mean(rr_for_init)
  # When the input is non-trivial (c_p or c_s != 0), the overall sample mean IBI
  # includes an input-driven mean-delta shift that biases the mu_0 anchor by
  # ≈ (c_s/a_s - c_p/a_p) * mean(u).  Use only near-baseline beats as anchor
  # when sufficient baseline data exist; otherwise fall back to the raw mean.
  mu_bar <- if (abs(params_init$free$c_p) + abs(params_init$free$c_s) > 1e-6) {
    u_at_beat_starts <- vapply(spikes[-length(spikes)], input_fn, numeric(1L))
    baseline_mask    <- abs(u_at_beat_starts) < 0.05
    if (sum(baseline_mask) >= 20L) mean(rr_for_init[baseline_mask])
    else mu_bar_raw
  } else {
    mu_bar_raw
  }

  params_spectral <- tryCatch(
    spectral_init(rr_for_init, verbose = verbose, mu_anchor = mu_bar),
    error = function(e) {
      if (verbose)
        message("pp_mle: spectral_init failed (", conditionMessage(e),
                "); falling back to params_init")
      params_init
    }
  )

  # Model structure (coupled vs uncoupled) is determined exclusively by params_init,
  # not by data-driven spectral estimates.  Band-filtered coupling initialisation is
  # applied only when params_init requests the coupled model.
  use_coupled <- (abs(params_init$free$a_ps) + abs(params_init$free$a_sp)) > 1e-12

  if (use_coupled) {
    cpl <- tryCatch(
      band_filtered_coupling_init(diff(spikes),
                                  a_p = params_spectral$free$a_p,
                                  a_s = params_spectral$free$a_s),
      error = function(e) list(a_ps = params_init$free$a_ps,
                               a_sp = params_init$free$a_sp,
                               c    = params_init$free$a_ps * params_init$free$a_sp)
    )

    params_spectral$free$a_ps <- max(cpl$a_ps, 1e-6)
    params_spectral$free$a_sp <- max(cpl$a_sp, 1e-6)
    # Stability guard: fall back if det(A) <= 0 at the starting point.
    if (params_spectral$free$a_p * params_spectral$free$a_s -
        params_spectral$free$a_ps * params_spectral$free$a_sp <= 0) {
      if (verbose)
        message("pp_mle: coupling starting values violate det(A)>0; using params_init.")
      params_spectral <- params_init
    }
  } else {
    params_spectral$free$a_ps <- 0
    params_spectral$free$a_sp <- 0
  }
  params_spectral$free$c_p <- params_init$free$c_p
  params_spectral$free$c_s <- params_init$free$c_s
  fp0 <- params_spectral$free

  # Optimizer coordinates: (log a_s, log(a_p − a_s), log A_p, log A_s, log κ [, log c])
  #
  # Coordinate | Expression         | Identified by
  # -----------|--------------------|----------------------------------------
  # log a_s    | log(a_s)           | LF ACF decay / spectral pole position
  # log gap    | log(a_p − a_s)     | HF pole position relative to LF pole
  # log A_p   | log(σ_p²/2a_p)    | Vagal spectral amplitude (OU stationary var of p)
  # log A_s   | log(σ_s²/2a_s)    | Sympathetic spectral amplitude (OU stationary var of s)
  # log κ     | log(μ₀/ρ²)        | IG shape; FIM diagonal in (log μ₀, log κ) — exact
  #
  # (log A_p, log A_s) replaces the previous (log σ_Δ², logit ψ) encoding.
  # Motivation (N. Beelders, pers. comm.): A_p = σ_p²/(2a_p) and A_s = σ_s²/(2a_s)
  # are the OU stationary variances and coincide with the Karhunen-Loève eigenvalue
  # envelopes of each branch.  When the poles are well-separated (a_p ≫ a_s,
  # the physiological regime), the two Lorentzians have negligible spectral overlap
  # and the Fisher information in (log A_p, log A_s) is approximately block-diagonal:
  # A_p is identified by HF variability, A_s by LF variability.  Equivalently, the
  # wavelet variance of the OU process at scale 1/a_ν is ≈ A_ν/2, so these
  # coordinates are the two dominant wavelet-variance coefficients of Δ(t).
  # The previous (log σ_Δ², logit ψ) encoding conflated the two amplitudes into
  # total + fraction, producing non-zero off-diagonal FIM terms that curved the
  # optimizer's descent path.
  #
  # Derived analytically in unpack():
  #   a_p   = a_s + exp(log_gap)          — a_p > a_s by construction
  #   σ_p   = √(2 a_p A_p),  σ_s = √(2 a_s A_s)   — exact, no approximation
  #   μ₀    = μ̄_obs · exp(−(A_p+A_s)/2)   — mean-anchor; exact under stationary model
  #
  # pack ∘ unpack is exact.  κ is the free coordinate; ρ is derived as sqrt(μ₀/κ).

  pack <- function(fp) {
    A_p  <- max(fp$sigma_p^2 / (2 * fp$a_p), 1e-12)
    A_s  <- max(fp$sigma_s^2 / (2 * fp$a_s), 1e-12)
    gap  <- max(fp$a_p - fp$a_s, 1e-4)
    mu_0_anchor <- max(mu_bar * exp(-(A_p + A_s) / 2), 0.10)
    kap_v <- kappa_from_rho(mu_0_anchor, fp$rho)   # = mu_0 / rho^2
    v <- c(log(fp$a_s),
           log(gap),
           log(A_p),
           log(A_s),
           log(max(kap_v, 0.5)))
    if (use_coupled) {
      c_val <- max(fp$a_ps * fp$a_sp, 1e-8)
      r_val <- max(fp$a_ps / max(fp$a_sp, 1e-8), 1e-3)
      v <- c(v, log(c_val), log(r_val))
    }
    v
  }

  unpack <- function(v) {
    a_s_v  <- exp(v[1L])
    a_p_v  <- a_s_v + exp(v[2L])        # a_p > a_s guaranteed
    A_p_v  <- max(exp(v[3L]), 1e-10)    # OU stationary variance of p
    A_s_v  <- max(exp(v[4L]), 1e-10)    # OU stationary variance of s
    kap_v  <- max(exp(v[5L]), 0.5)
    mu_0_v <- max(mu_bar * exp(-(A_p_v + A_s_v) / 2), 0.10)
    rho_v  <- sqrt(mu_0_v / kap_v)   # derived, not free
    c_v   <- if (use_coupled) exp(v[6L]) else 0
    r_v   <- if (use_coupled) exp(v[7L]) else 1
    a_ps_v <- sqrt(max(c_v * r_v, 0))
    a_sp_v <- sqrt(max(c_v / r_v, 0))
    list(
      a_p     = a_p_v,
      a_s     = a_s_v,
      a_ps    = a_ps_v,
      a_sp    = a_sp_v,
      sigma_p = sqrt(max(2 * a_p_v * A_p_v, 1e-10)),
      sigma_s = sqrt(max(2 * a_s_v * A_s_v, 1e-10)),
      mu_0    = mu_0_v,
      rho     = rho_v,
      c_p     = fp0$c_p,
      c_s     = fp0$c_s
    )
  }

  # Bounds in optimizer coordinates.
  # v[1] log(a_s):   a_s ∈ [0.01, 2.0] Hz
  # v[2] log(gap):   gap ∈ [0.01, 7.99] → a_p ≤ 9.99 Hz
  # v[3] log(A_p):   A_p ∈ [1e-6, 0.45]  (vagal stationary variance)
  # v[4] log(A_s):   A_s ∈ [1e-6, 0.45]  (sympathetic stationary variance)
  # v[5] log(kappa):     kappa ∈ [0.5, 1000]
  # Individual caps at 0.45 keep total σ_Δ² = A_p + A_s well below 1,
  # matching the physiological constraint that OU variance ≪ μ₀².
  lower_v <- c(log(0.01), log(0.01), log(1e-6), log(1e-6), log(0.5))
  upper_v <- c(log(2.0),  log(7.99), log(0.45), log(0.45), log(1000))

  if (use_coupled) {
    lower_v <- c(lower_v, log(1e-8),  log(1e-3))   # log(c), log(r)
    # No fixed upper bound on log(c): the stability constraint a_p*a_s - a_ps*a_sp > 0
    # is enforced analytically inside neg_ll for every candidate parameter vector.
    # A static bound derived from the initial (fp0$a_p, fp0$a_s) is incorrect once
    # the optimizer moves away from those starting values.
    upper_v <- c(upper_v,
                 log(0.99 * max(fp0$a_p, 2.0) * max(fp0$a_s, 0.5)),  # generous ceiling
                 log(1e3))
  }

  neg_ll <- function(v) {
    fp_v <- unpack(v)
    # a_p > a_s is guaranteed by construction: a_p = a_s + exp(log_gap).
    # For the coupled model, additionally enforce det(A) = a_p*a_s − a_ps*a_sp > 0.
    if (use_coupled && fp_v$a_p * fp_v$a_s - fp_v$a_ps * fp_v$a_sp <= 0) {
      return(1e10)
    }

    params <- list(structural = params_init$structural, free = fp_v)
    flt    <- tryCatch(pp_ukf(spikes, params, input_fn),
                       error = function(e) list(innov = NULL))
    if (is.null(flt$innov)) return(1e10)
    ll_vec_fin <- flt$ll_vec[is.finite(flt$ll_vec)]
    if (length(ll_vec_fin) < 2L) return(1e10)
    -sum(ll_vec_fin)
  }

  res <- optim(
    par     = pack(fp0),
    fn      = neg_ll,
    method  = "L-BFGS-B",
    lower   = lower_v,
    upper   = upper_v,
    control = list(maxit = 10000L,
                   factr = 1e6,
                   trace = if (verbose) 1L else 0L)
  )

  if (res$convergence != 0L && verbose)
    message("pp_mle: L-BFGS-B convergence code ", res$convergence,
            " — consider re-running from a perturbed starting value.")

  fp_hat     <- unpack(res$par)
  params_hat <- list(structural = params_init$structural, free = fp_hat)

  ll_hat <- if (res$value >= 1e9) {
    if (verbose)
      message("pp_mle: optimizer terminated at a boundary (neg_ll = ", res$value,
              "). ll is unreliable; check convergence code ", res$convergence, ".")
    NA_real_
  } else {
    -res$value
  }

  list(
    params_hat  = params_hat,
    free_hat    = fp_hat,
    ll          = ll_hat,
    convergence = res$convergence,
    filter      = pp_ukf(spikes, params_hat, input_fn)
  )
}

# ---- Concentrated MLE: poles fixed, amplitudes and shape free ----
#
# Reduces the 5D pp_mle optimisation to 3D by fixing (a_p, a_s) at the
# supplied values (typically from wavelet_spectral_fit Stage 1).
# Free coordinates: (log A_p, log A_s, log κ).
# Mean anchor μ₀ = μ̄·exp(−(A_p+A_s)/2) is retained exactly as in pp_mle.
# Returns the same structure as pp_mle() for drop-in compatibility.

pp_mle_concentrated <- function(spikes,
                                params_init,
                                a_p_fixed,
                                a_s_fixed,
                                input_fn = function(t) 0,
                                verbose  = FALSE) {
  rr <- diff(spikes)
  if (length(rr) < 10L) stop("pp_mle_concentrated: need at least 10 spikes.")
  mu_bar_c_raw <- mean(rr)
  mu_bar_c <- if (abs(params_init$free$c_p) + abs(params_init$free$c_s) > 1e-6) {
    u_at_beat_starts <- vapply(spikes[-length(spikes)], input_fn, numeric(1L))
    baseline_mask    <- abs(u_at_beat_starts) < 0.05
    if (sum(baseline_mask) >= 20L) mean(rr[baseline_mask]) else mu_bar_c_raw
  } else {
    mu_bar_c_raw
  }

  # Starting values: spectral_init re-evaluated at the fixed poles
  sp_init <- tryCatch(spectral_init(rr, verbose = FALSE),
                      error = function(e) params_init)
  fp_sp   <- sp_init$free
  A_p0    <- max(fp_sp$sigma_p^2 / (2 * a_p_fixed), 1e-8)
  A_s0    <- max(fp_sp$sigma_s^2 / (2 * a_s_fixed), 1e-8)
  kap0    <- kappa_from_rho(fp_sp$mu_0, fp_sp$rho)
  v0      <- c(log(max(A_p0, 1e-8)), log(max(A_s0, 1e-8)), log(max(kap0, 0.5)))

  neg_ll_c <- function(v) {
    A_p_v  <- max(exp(v[1L]), 1e-10)
    A_s_v  <- max(exp(v[2L]), 1e-10)
    kap_v  <- max(exp(v[3L]), 0.5)
    mu_0_v <- max(mu_bar_c * exp(-(A_p_v + A_s_v) / 2), 0.10)
    fp_v <- list(
      a_p     = a_p_fixed,
      a_s     = a_s_fixed,
      a_ps    = params_init$free$a_ps,
      a_sp    = params_init$free$a_sp,
      sigma_p = sqrt(max(2 * a_p_fixed * A_p_v, 1e-10)),
      sigma_s = sqrt(max(2 * a_s_fixed * A_s_v, 1e-10)),
      mu_0    = mu_0_v,
      rho     = sqrt(mu_0_v / kap_v),
      c_p     = params_init$free$c_p,
      c_s     = params_init$free$c_s
    )
    params_v <- list(structural = params_init$structural, free = fp_v)
    flt <- tryCatch(pp_ukf(spikes, params_v, input_fn),
                    error = function(e) list(innov = NULL))
    if (is.null(flt$innov)) return(1e10)
    ll_vec_fin <- flt$ll_vec[is.finite(flt$ll_vec)]
    if (length(ll_vec_fin) < 2L) return(1e10)
    -sum(ll_vec_fin)
  }

  res <- optim(v0, neg_ll_c,
               method  = "L-BFGS-B",
               lower   = c(log(1e-6), log(1e-6), log(0.5)),
               upper   = c(log(0.45), log(0.45), log(1000)),
               control = list(maxit = 10000, factr = 1e6,
                              trace = if (verbose) 1L else 0L))

  if (res$convergence != 0L && verbose)
    message("pp_mle_concentrated: L-BFGS-B code ", res$convergence)

  th     <- res$par
  A_p_h  <- exp(th[1L]);  A_s_h <- exp(th[2L]);  kap_h <- exp(th[3L])
  mu_0_h <- max(mu_bar_c * exp(-(A_p_h + A_s_h) / 2), 0.10)

  fp_hat <- list(
    a_p     = a_p_fixed,
    a_s     = a_s_fixed,
    a_ps    = params_init$free$a_ps,
    a_sp    = params_init$free$a_sp,
    sigma_p = sqrt(max(2 * a_p_fixed * A_p_h, 1e-10)),
    sigma_s = sqrt(max(2 * a_s_fixed * A_s_h, 1e-10)),
    mu_0    = mu_0_h,
    rho     = sqrt(mu_0_h / kap_h),
    c_p     = params_init$free$c_p,
    c_s     = params_init$free$c_s
  )
  params_hat <- list(structural = params_init$structural, free = fp_hat)
  ll_hat     <- if (res$value >= 1e9) NA_real_ else -res$value

  list(params_hat  = params_hat,
       free_hat    = fp_hat,
       ll          = ll_hat,
       convergence = res$convergence,
       filter      = pp_ukf(spikes, params_hat, input_fn))
}

# ---- Two-stage concentrated MLE ----
#
# Stage 1 — wavelet_spectral_fit() [spectral_init.R]:
#   Estimates (a_p, a_s) from tachogram Haar level variances via the Whittle
#   likelihood. Globally convex; always finds the global optimum.
#
# Stage 2 — pp_mle_concentrated():
#   Fixes Stage-1 poles; optimises 3D (log A_p, log A_s, log κ).
#   Eliminates the σ–a ridge and μ₀–σ_Δ² cross-block coupling.
#
# Stage 3 (optional, refine = TRUE) — pp_mle():
#   Full 5D refinement from the Stage-2 solution; keeps the better LL.
#   Set refine = FALSE for speed-critical loops (sensitivity, large-N recovery).
#
# Falls back to pp_mle() at any stage failure.
# Returns the same structure as pp_mle() for drop-in compatibility.

pp_mle_twostage <- function(spikes,
                            params_init,
                            input_fn = function(t) 0,
                            refine   = TRUE,
                            verbose  = FALSE) {
  rr <- diff(spikes)

  # Stage 1
  wave <- tryCatch(
    wavelet_spectral_fit(rr),
    error = function(e) {
      if (verbose) message("pp_mle_twostage: Stage 1 failed (", conditionMessage(e),
                           "); using pp_mle.")
      NULL
    }
  )
  if (is.null(wave) || !is.finite(wave$a_p) || !is.finite(wave$a_s)) {
    if (verbose) message("pp_mle_twostage: Stage 1 invalid; using pp_mle.")
    return(pp_mle(spikes, params_init, input_fn, verbose = verbose))
  }
  if (verbose) message(sprintf(
    "pp_mle_twostage: Stage 1  a_p=%.3f  a_s=%.4f  A_p=%.4f  A_s=%.4f",
    wave$a_p, wave$a_s, wave$A_p, wave$A_s))

  # Stage 2
  conc <- tryCatch(
    pp_mle_concentrated(spikes, params_init,
                        a_p_fixed = wave$a_p, a_s_fixed = wave$a_s,
                        input_fn  = input_fn,  verbose   = verbose),
    error = function(e) {
      if (verbose) message("pp_mle_twostage: Stage 2 failed (", conditionMessage(e),
                           "); using pp_mle.")
      NULL
    }
  )
  if (is.null(conc))
    return(pp_mle(spikes, params_init, input_fn, verbose = verbose))

  if (verbose && is.finite(conc$ll))
    message(sprintf("pp_mle_twostage: Stage 2  ll=%.2f  (convergence=%d)",
                    conc$ll, conc$convergence))

  # Stage 3 (optional)
  if (!refine) return(conc)

  full <- tryCatch(
    pp_mle(spikes, conc$params_hat, input_fn, verbose = verbose),
    error = function(e) NULL
  )

  if (is.null(full) || is.na(full$ll) ||
      (!is.na(conc$ll) && conc$ll >= full$ll)) {
    if (verbose && !is.null(full) && is.finite(full$ll))
      message(sprintf("pp_mle_twostage: concentrated retained (%.2f >= %.2f)",
                      conc$ll, full$ll))
    return(conc)
  }
  full
}

# ---- Interpolate filtered state to a regular time grid ----
# Useful for overlaying filter output on continuous simulation plots.

# OU-propagated inter-beat interpolation helper
.ou_propagate_grid <- function(beat_times, m_beat, a_val, time_grid) {
  # findInterval(x, v) returns k s.t. v[k] <= x < v[k+1], 0 before first.
  grp   <- findInterval(time_grid, beat_times)
  n_bt  <- length(beat_times)
  valid <- grp >= 1L & grp <= n_bt
  out   <- rep(NA_real_, length(time_grid))
  if (!any(valid)) return(out)
  gv         <- grp[valid]
  out[valid] <- m_beat[gv] * exp(-a_val * (time_grid[valid] - beat_times[gv]))
  out
}

.prop_sd_delta <- function(ukf_result, time_grid, fp) {
  bt  <- ukf_result$beat_times
  out <- rep(NA_real_, length(time_grid))
  n_bt <- length(bt)

  for (k in seq_len(n_bt)) {
    mask <- if (k < n_bt) time_grid >= bt[k] & time_grid < bt[k + 1L]
    else           time_grid >= bt[k]
    if (!any(mask)) next

    dt_v <- time_grid[mask] - bt[k]
    P_k  <- ukf_result$P_filt[[k]]

    F_p <- exp(-fp$a_p * dt_v)
    F_s <- exp(-fp$a_s * dt_v)
    Q_p <- (fp$sigma_p^2 / (2 * fp$a_p)) * (-expm1(-2 * fp$a_p * dt_v))
    Q_s <- (fp$sigma_s^2 / (2 * fp$a_s)) * (-expm1(-2 * fp$a_s * dt_v))

    var_p  <- F_p^2 * P_k[1L, 1L] + Q_p
    var_s  <- F_s^2 * P_k[2L, 2L] + Q_s
    cov_ps <- F_p * F_s * P_k[1L, 2L]

    out[mask] <- sqrt(pmax(var_p + var_s - 2 * cov_ps, 0))
  }
  out
}

filter_to_grid <- function(ukf_result, time_grid) {
  bt   <- ukf_result$beat_times
  fp   <- ukf_result$params$free

  p_prop <- .ou_propagate_grid(bt, ukf_result$m_filt[, "p"],
                               fp$a_p, time_grid)
  s_prop <- .ou_propagate_grid(bt, ukf_result$m_filt[, "s"],
                               fp$a_s, time_grid)
  first_bt <- bt[1L]
  p_prop[time_grid < first_bt] <- NA_real_
  s_prop[time_grid < first_bt] <- NA_real_

  interp <- function(y) {              # keep linear for scalar quantities
    out <- approx(bt, y, xout = time_grid, rule = 2)$y
    out[time_grid < first_bt] <- NA_real_
    out
  }

  list(
    delta    = s_prop - p_prop,
    p        = p_prop,
    s        = s_prop,
    mu       = interp(ukf_result$mu_filt),
    sd_delta = .prop_sd_delta(ukf_result, time_grid, fp)
  )
}

# Helper — Propagate per-branch posterior variance with OU dynamics
.prop_sd_branch <- function(ukf_result, time_grid, fp, branch) {
  # branch = 1 (p) or 2 (s)
  bt   <- ukf_result$beat_times
  a_v  <- if (branch == 1L) fp$a_p else fp$a_s
  s2_v <- if (branch == 1L) fp$sigma_p^2 else fp$sigma_s^2
  out  <- rep(NA_real_, length(time_grid))
  for (k in seq_len(length(bt))) {
    mask <- if (k < length(bt)) time_grid >= bt[k] & time_grid < bt[k + 1L]
    else                time_grid >= bt[k]
    if (!any(mask)) next
    dt_v  <- time_grid[mask] - bt[k]
    P_kbb <- ukf_result$P_filt[[k]][branch, branch]
    Q_v   <- (s2_v / (2 * a_v)) * (-expm1(-2 * a_v * dt_v))
    out[mask] <- sqrt(pmax(exp(-2 * a_v * dt_v) * P_kbb + Q_v, 0))
  }
  out
}
