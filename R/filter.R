# ============================================================
# R/filter.R
#
# Spike-to-spike Unscented Kalman Filter for the SDE-IG model.
# Requires: model_functions.R sourced first.
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
  L <- 2L; alpha <- 0.1; beta <- 2.0; kappa <- 0.0
  lambda  <- alpha^2 * (L + kappa) - L
  c_val   <- L + lambda
  n_pts   <- 2L * L + 1L
  Wm <- c(lambda / c_val, rep(0.5 / c_val, 2L * L))
  Wc <- c(lambda / c_val + (1 - alpha^2 + beta), rep(0.5 / c_val, 2L * L))
  list(lambda = lambda, c = c_val, Wm = Wm, Wc = Wc, L = L, n_pts = n_pts)
})

# ---- Sigma-point generation (scaled, Cholesky-based) ----

.sigma_points <- function(m, P) {
  L <- .ukf_weights$L
  A <- tryCatch(
    t(chol(.ukf_weights$c * P)),
    error = function(e)
      t(chol(.ukf_weights$c * (P + 1e-8 * diag(L))))  # regularize if needed
  )
  pts <- matrix(0, L, .ukf_weights$n_pts)
  pts[, 1L] <- m
  for (i in seq_len(L)) {
    pts[, i + 1L]      <- m + A[, i]
    pts[, i + 1L + L]  <- m - A[, i]
  }
  pts
}

# ---- Exact OU prediction ----
# The OU transition is linear-Gaussian, so the prediction step is an
# exact Kalman prediction — no approximation.
#
# Returns: list(m, P)

.ou_predict <- function(m0, P0, tau, sp, fp, u = 0) {
  a_p <- sp$a_p;  a_s <- sp$a_s
  F_p <- exp(-a_p * tau);  F_s <- exp(-a_s * tau)

  # Input-driven mean shift (exact OU integral of the forcing term)
  d_p <- (fp$c_p * u / a_p) * (-expm1(-a_p * tau))
  d_s <- (fp$c_s * u / a_s) * (-expm1(-a_s * tau))

  m1 <- c(m0[1L] * F_p + d_p,
          m0[2L] * F_s + d_s)

  # Covariance: F P Fᵀ + Q  (F diagonal ⇒ off-diagonal of F P Fᵀ = F_p F_s P[1,2])
  Q_p <- (fp$sigma_p^2 / (2 * a_p)) * (-expm1(-2 * a_p * tau))
  Q_s <- (fp$sigma_s^2 / (2 * a_s)) * (-expm1(-2 * a_s * tau))

  P1       <- matrix(0, 2, 2)
  P1[1, 1] <- F_p^2 * P0[1, 1] + Q_p
  P1[2, 2] <- F_s^2 * P0[2, 2] + Q_s
  P1[1, 2] <- F_p * F_s * P0[1, 2]
  P1[2, 1] <- P1[1, 2]

  list(m = m1, P = P1)
}

# ---- UKF update at a spike ----
# Observation: τ_k ~ IG(μ₀ exp(p − s), κ).
# Observation function: h(x) = μ₀ exp(p − s).
# Noise variance: R = μ̂³/κ (IG variance outer-linearized at predicted mean).
#
# Returns: list(m, P, innov, S, mu_hat, ll)

.ukf_update <- function(m1, P1, tau_k, fp) {
  mu_0  <- fp$mu_0
  kappa <- fp$kappa

  pts    <- .sigma_points(m1, P1)
  h_pts  <- mu_0 * exp(pts[1L, ] - pts[2L, ])   # h(sigma_i)

  mu_hat <- sum(.ukf_weights$Wm * h_pts)
  R      <- max(mu_hat^3 / kappa, 1e-6)          # IG Var; floored for safety

  dh   <- h_pts - mu_hat
  S    <- sum(.ukf_weights$Wc * dh^2) + R

  # Cross-covariance state–observation (2×1)
  dx   <- pts - m1                                 # 2 × n_pts
  C    <- rowSums(sweep(dx, 2L, .ukf_weights$Wc * dh, `*`))

  K      <- C / S
  innov  <- tau_k - mu_hat

  m2     <- m1 + K * innov
  P2_raw <- P1 - outer(K, K) * S
  P2     <- 0.5 * (P2_raw + t(P2_raw)) + 1e-9 * diag(.ukf_weights$L)

  # Marginal LL: Gaussian innovation approximation (standard UKF)
  ll_k <- dnorm(tau_k, mean = mu_hat, sd = sqrt(max(S, 1e-9)), log = TRUE)

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
  N  <- length(spikes)
  if (N < 2L) stop("pp_ukf: need at least 2 spikes.")

  if (is.null(m0)) m0 <- c(0, 0)
  if (is.null(P0)) P0 <- diag(c(fp$sigma_p^2 / (2 * sp$a_p),
                                fp$sigma_s^2 / (2 * sp$a_s)))

  n_ibi      <- N - 1L
  m_filt     <- matrix(0, n_ibi, 2L, dimnames = list(NULL, c("p", "s")))
  P_filt     <- vector("list", n_ibi)
  delta_filt <- mu_filt <- innov_v <- S_v <- ll_v <- numeric(n_ibi)

  m_cur <- m0; P_cur <- P0

  for (k in seq_len(n_ibi)) {
    tau_k <- spikes[k + 1L] - spikes[k]
    u_k   <- input_fn(spikes[k])     # constant-input approximation over IBI

    pred <- .ou_predict(m_cur, P_cur, tau_k, sp, fp, u = u_k)
    upd  <- .ukf_update(pred$m, pred$P, tau_k, fp)

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
                   input_fn       = function(t) 0,
                   optimize_gains = FALSE,
                   verbose        = FALSE) {
  sp <- params_init$structural

  # Log-space packing: (log σ_p, log σ_s, log μ₀, log κ [, c_p, c_s])
  pack <- function(fp) {
    v <- c(log(fp$sigma_p), log(fp$sigma_s), log(fp$mu_0), log(fp$kappa))
    if (optimize_gains) c(v, fp$c_p, fp$c_s) else v
  }
  unpack <- function(v) {
    list(sigma_p = exp(v[1L]),
         sigma_s = exp(v[2L]),
         mu_0    = exp(v[3L]),
         kappa   = exp(v[4L]),
         c_p     = if (optimize_gains) v[5L] else params_init$free$c_p,
         c_s     = if (optimize_gains) v[6L] else params_init$free$c_s)
  }

  neg_ll <- function(v) {
    fp_v   <- unpack(v)
    params <- list(structural = sp, free = fp_v)
    flt    <- tryCatch(pp_ukf(spikes, params, input_fn),
                       error = function(e) list(ll = -Inf))
    if (!is.finite(flt$ll)) return(1e10)
    -flt$ll
  }

  res <- optim(
    par     = pack(params_init$free),
    fn      = neg_ll,
    method  = "L-BFGS-B",
    control = list(maxit = 500L,
                   factr = 1e7,            # moderate tolerance; tighten for final run
                   trace = if (verbose) 1L else 0L)
  )

  if (res$convergence != 0L && verbose)
    message("pp_mle: L-BFGS-B convergence code ", res$convergence,
            " — consider re-running from a perturbed starting value.")

  fp_hat     <- unpack(res$par)
  params_hat <- list(structural = sp, free = fp_hat)

  list(
    params_hat  = params_hat,
    free_hat    = fp_hat,
    ll          = -res$value,
    convergence = res$convergence,
    filter      = pp_ukf(spikes, params_hat, input_fn)
  )
}

# ---- Interpolate filtered state to a regular time grid ----
# Useful for overlaying filter output on continuous simulation plots.

filter_to_grid <- function(ukf_result, time_grid) {
  bt <- ukf_result$beat_times
  list(
    delta = approx(bt, ukf_result$delta_filt,          xout = time_grid, rule = 2)$y,
    p     = approx(bt, ukf_result$m_filt[, "p"],       xout = time_grid, rule = 2)$y,
    s     = approx(bt, ukf_result$m_filt[, "s"],       xout = time_grid, rule = 2)$y,
    mu    = approx(bt, ukf_result$mu_filt,              xout = time_grid, rule = 2)$y,
    sd_delta = approx(bt,
                      sqrt(sapply(ukf_result$P_filt,
                                  function(P) P[1,1] + P[2,2] - 2*P[1,2])),
                      xout = time_grid, rule = 2)$y
  )
}
