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

  Q_p <- (fp$sigma_p^2 / (2 * a_p)) * (1 - F_p^2)   # exact: 1-e^{-2at} = (1-F)(1+F)
  Q_s <- (fp$sigma_s^2 / (2 * a_s)) * (1 - F_s^2)
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
    # For the double-logistic protocol (10--90% rise ≈ 0.44 s) and IBIs ~0.85 s,
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
  params_spectral <- tryCatch(
    spectral_init(rr_for_init, verbose = verbose),
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
                               a_sp = params_init$free$a_sp)
    )
    # Clamp away from the L-BFGS-B lower boundary so gradients are finite at init.
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

  pack <- function(fp) {
    v <- c(log(fp$a_p), log(fp$a_s),
           log(fp$sigma_p), log(fp$sigma_s),
           log(fp$mu_0), log(fp$rho))
    if (use_coupled) v <- c(v, fp$a_ps, fp$a_sp)
    v
  }
  unpack <- function(v) {
    list(
      a_p     = exp(v[1L]),
      a_s     = exp(v[2L]),
      a_ps    = if (use_coupled) max(v[7L], 0) else 0,
      a_sp    = if (use_coupled) max(v[8L], 0) else 0,
      sigma_p = exp(v[3L]),
      sigma_s = exp(v[4L]),
      mu_0    = exp(v[5L]),
      rho     = exp(v[6L]),
      c_p     = fp0$c_p,  # fixed known constant carried from params_init
      c_s     = fp0$c_s   # not estimated; input contribution analytically removed
    )
  }

  # Bounds reflect the physiological range of each parameter.
  # a_p in [0.30, 10.0] Hz  (vagal time constant 0.1 – 3 s)
  # a_s in [0.01,  2.0] Hz  (sympathetic time constant 0.5 – 100 s)
  # Tight bounds prevent the optimizer from sliding along the marginal
  # identifiability ridge  sigma_p^2/(2*a_p) = const  to implausible
  # values such as a_p ~ 50 Hz.
  lower_v <- c(log(0.30), log(0.01), log(1e-4), log(1e-4), log(0.25), log(0.01))
  upper_v <- c(log(10.0), log(2.0),  log(10.0), log(10.0), log(4.0),  log(2.0))

  if (use_coupled) {
    lower_v <- c(lower_v, 0, 0)    # a_ps, a_sp >= 0
    upper_v <- c(upper_v, 10, 10)  # coupling terms <= 10 (stability enforced separately)
  }

  neg_ll <- function(v) {
    fp_v <- unpack(v)
    # Hard constraint: a_p > a_s ensures LF/HF pole separation and prevents
    # the optimizer from swapping branch identities.
    if (fp_v$a_p <= fp_v$a_s) return(1e10)
    # Coupled-OU stability: det(A) > 0.
    if (fp_v$a_p * fp_v$a_s - fp_v$a_ps * fp_v$a_sp <= 0) return(1e10)
    params <- list(structural = params_init$structural, free = fp_v)
    flt    <- tryCatch(pp_ukf(spikes, params, input_fn),
                       error = function(e) list(innov = NULL))
    if (is.null(flt$innov)) return(1e10)
    # Gaussian innovation log-likelihood  log N(innov; 0, S_k).
    # Provides continuous, well-defined gradients for L-BFGS-B throughout
    # the parameter space.  The sigma-point IG average stored in flt$ll is
    # reserved for state-estimation diagnostics; it can exhibit near-
    # discontinuous Jacobians when sigma-points straddle the IG hazard peak,
    # degrading gradient-based optimization.
    fin <- is.finite(flt$innov) & is.finite(flt$S_innov) & flt$S_innov > 0
    if (!any(fin)) return(1e10)
    ll_g <- sum(dnorm(flt$innov[fin], 0, sqrt(flt$S_innov[fin]), log = TRUE))
    if (!is.finite(ll_g)) return(1e10)
    -ll_g
  }

  res <- optim(
    par     = pack(fp0),
    fn      = neg_ll,
    method  = "L-BFGS-B",
    lower   = lower_v,
    upper   = upper_v,
    control = list(maxit = 1000L,
                   factr = 1e6,     # tighter convergence tolerance
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

# ---- Interpolate filtered state to a regular time grid ----
# Useful for overlaying filter output on continuous simulation plots.

filter_to_grid <- function(ukf_result, time_grid) {
  bt <- ukf_result$beat_times

  first_bt <- bt[1L]

  interp <- function(y) {
    out <- approx(bt, y, xout = time_grid, rule = 2)$y
    out[time_grid < first_bt] <- NA_real_
    out
  }

  list(
    delta    = interp(ukf_result$delta_filt),
    p        = interp(ukf_result$m_filt[, "p"]),
    s        = interp(ukf_result$m_filt[, "s"]),
    mu       = interp(ukf_result$mu_filt),
    sd_delta = interp(sqrt(pmax(vapply(ukf_result$P_filt,
                                       function(P) P[1L,1L] + P[2L,2L] - 2*P[1L,2L],
                                       numeric(1L)), 0)))
  )
}
