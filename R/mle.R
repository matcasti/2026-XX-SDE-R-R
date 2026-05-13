# ============================================================
# R/mle.R
#
# Closed-form and numerical MLE for the SDE-IG model.
# Requires: model_functions.R sourced first.
#
# KEY RESULT: the conditional log-likelihood (given latent states)
# separates into three independent blocks:
#
#   L(theta | states) = L_p(sigma_p, c_p)    [OU parasympathetic]
#                     + L_s(sigma_s, c_s)    [OU sympathetic]
#                     + L_obs(mu_0, kappa)   [IG observations]
#
# Each block has a closed-form MLE. The FIM is therefore block-diagonal,
# and positive-definiteness of each block proves identifiability.
# ============================================================


# Helper to compute the effective integrate-and-fire drift
compute_effective_delta <- function(spk, time_vec, delta_vec) {
  # findInterval gives O(log n_time) per IBI vs O(n_time) from which().
  # time_vec is strictly sorted and uniform; semantics match the original:
  #   i_start = last index with time_vec[i] <= t_start
  #   i_end   = last index with time_vec[i] <= t_end
  #   qualifying indices: (i_start+1):i_end  ↔  time_vec > t_start & <= t_end
  n_ibi    <- length(spk) - 1L
  i_starts <- findInterval(spk[-length(spk)], time_vec)
  i_ends   <- findInterval(spk[-1L],          time_vec)
  dt_val   <- if (length(time_vec) > 1L) time_vec[2L] - time_vec[1L] else NA_real_

  vapply(seq_len(n_ibi), function(k) {
    i0 <- i_starts[k];  i1 <- i_ends[k]
    if (i1 <= i0) {
      warning(sprintf(
        "compute_effective_delta: IBI [%.4f, %.4f] (%.4f s) < dt (%.4f s); returning NA.",
        spk[k], spk[k + 1L], spk[k + 1L] - spk[k], dt_val))
      return(NA_real_)
    }
    idx   <- (i0 + 1L):i1
    neg_d <- -delta_vec[idx]
    d_max <- max(neg_d)
    log(length(idx)) - (d_max + log(sum(exp(neg_d - d_max))))
  }, numeric(1L))
}

# ---- Block 1 / 2: Analytic OU MLE ----
#
# Model: x[i+1] | x[i] ~ N(x[i]*exp(-a*dt) + (c*u[i]/a)*(1-exp(-a*dt)),
#                           sigma^2/(2a)*(1-exp(-2a*dt)))
#
# Closed-form solution from OLS (for c) + method-of-moments (for sigma).

ou_mle <- function(x_vec, u_vec, a_init, dt,
                   log_a_bounds = c(log(0.01), log(100)),
                   c_fixed = NULL) {
  n <- length(x_vec) - 1L
  if (n < 2L)
    return(list(a = a_init, sigma = NA_real_, c_gain = NA_real_,
                a_se = NA_real_, sigma_se = NA_real_, c_se = NA_real_,
                log_lik = NA_real_, n_obs = n))

  idx    <- seq_len(n)
  x_cur  <- x_vec[idx]          # x_1 … x_{n}   — stable across all profile_at_a calls
  x_next <- x_vec[-1L]          # x_2 … x_{n+1}
  u_cur  <- u_vec[idx]          # u_1 … u_{n}

  profile_at_a <- function(a) {
    e_adt  <- exp(-a * dt)
    z      <- u_cur / a * (-expm1(-a * dt))
    y      <- x_next - x_cur * e_adt
    C_a    <- (-expm1(-2 * a * dt)) / (2 * a)
    c_hat  <- if (!is.null(c_fixed)) {
      c_fixed
    } else {
      ssz <- sum(z^2)
      if (ssz > 1e-14) {
        sum(y * z) / ssz
      } else { 0 }
    }
    resid  <- y - c_hat * z
    sigma2 <- max(mean(resid^2) / C_a, 1e-14)
    ll     <- sum(dnorm(resid, 0, sqrt(sigma2 * C_a), log = TRUE))
    list(a = a, c_gain = c_hat, sigma = sqrt(sigma2),
         ll = ll, C_a = C_a)
  }

  # 1-D profile optimisation over log(a)
  last_pf <- NULL
  opt <- tryCatch(
    optimize(function(la) {
      pf <- profile_at_a(exp(la))
      if (is.finite(pf$ll)) { last_pf <<- pf; pf$ll } else -Inf
    }, interval = log_a_bounds, maximum = TRUE, tol = 1e-7),
    error = function(e) list(maximum = log(a_init), objective = -Inf))
  a_hat <- exp(opt$maximum)
  pf    <- if (!is.null(last_pf)) last_pf else profile_at_a(a_hat)

  # Numerical Hessian for standard errors.
  # When u(t) ≡ 0, nll_3d does not depend on th[3] (c_gain is non-identifiable),
  # making the 3×3 Hessian singular. Detect this case and use a 2×2 Hessian
  # in (log a, log σ) only; c_se is set to NA (not estimable).
  th0    <- c(log(a_hat), log(pf$sigma), pf$c_gain)

  # c is not a free parameter when c_fixed is supplied OR when u ≡ 0.
  u_is_zero <- !is.null(c_fixed) || all(abs(u_cur) < .Machine$double.eps * 100)

  nll_3d <- function(th) {
    a_v   <- exp(th[1L]); sig_v <- exp(th[2L]); c_v <- th[3L]
    if (a_v < 1e-6 || sig_v <= 0) return(1e10)
    e_adt <- exp(-a_v * dt)
    z_v   <- u_cur / a_v * (-expm1(-a_v * dt))        # use pre-indexed u_cur
    y_v   <- x_next - x_cur * e_adt - c_v * z_v       # use pre-indexed x_cur, x_next
    C_v   <- (-expm1(-2 * a_v * dt)) / (2 * a_v)
    ll_v  <- tryCatch(sum(dnorm(y_v, 0, sig_v * sqrt(C_v), log = TRUE)),
                      error = function(e) -Inf)
    if (is.finite(ll_v)) -ll_v else 1e10
  }

  if (u_is_zero) {
    # 2×2 Hessian in (log a, log σ); c is absorbed into the residual variance.
    th0_2 <- th0[1:2]
    nll_2 <- function(th) nll_3d(c(th, if (!is.null(c_fixed)) c_fixed else 0))
    H2      <- numerical_hessian(nll_2, th0_2)
    V2      <- tryCatch(solve(H2), error = function(e) matrix(NA_real_, 2L, 2L))
    safe2   <- function(k) sqrt(max(V2[k, k], 0, na.rm = TRUE))
    return(list(
      a        = a_hat,
      sigma    = pf$sigma,
      c_gain   = if (!is.null(c_fixed)) c_fixed else 0,
      a_se     = a_hat    * safe2(1L),
      sigma_se = pf$sigma * safe2(2L),
      c_se     = NA_real_,
      log_lik  = pf$ll,
      n_obs    = n
    ))
  }

  H <- numerical_hessian(nll_3d, th0)

  V        <- tryCatch(solve(H), error = function(e) matrix(NA_real_, 3L, 3L))
  safe_se  <- function(k) sqrt(max(V[k, k], 0, na.rm = TRUE))

  list(
    a        = a_hat,
    sigma    = pf$sigma,
    c_gain   = pf$c_gain,
    a_se     = a_hat     * safe_se(1L),
    sigma_se = pf$sigma  * safe_se(2L),
    c_se     = safe_se(3L),
    log_lik  = pf$ll,
    n_obs    = n
  )
}

# ---- Bivariate coupled-OU log-likelihood ----
# Parameterization: all four rate parameters free (a_p, a_s, a_ps, a_sp);
# called by ou_coupled_mle which optimises over log(a_p), log(a_s), a_ps, a_sp.
# All n-1 bivariate Gaussian transition densities share the same F, Q
# (dt constant), so these are computed once per call.

ou_coupled_log_lik <- function(a_ps, a_sp, a_p, a_s, sigma_p, sigma_s, c_p, c_s,
                               x_mat, u_vec, dt) {
  if (a_p * a_s - a_ps * a_sp <= 0) return(-Inf)
  mats <- tryCatch(
    ou_coupled_matrices(a_p, a_s, a_ps, a_sp, sigma_p, sigma_s, dt),
    error = function(e) NULL
  )
  if (is.null(mats) || !is.finite(mats$log_det_Q)) return(-Inf)

  n      <- nrow(x_mat) - 1L
  c_vec  <- c(c_p, c_s)
  ll_base <- -n * log(2 * pi) - (n / 2) * mats$log_det_Q

  X_cur  <- x_mat[-nrow(x_mat), , drop = FALSE]   # n × 2
  X_next <- x_mat[-1L,          , drop = FALSE]   # n × 2

  means  <- X_cur %*% mats$Ft

  if (any(u_vec[seq_len(n)] != 0)) {
    d_base <- as.vector(mats$FmI %*% mats$A_inv %*% c_vec)
    means  <- means + outer(u_vec[seq_len(n)], d_base)
  }

  resids <- X_next - means                         # n × 2
  # Hadamard trace trick: sum_i rᵢᵀ Q⁻¹ rᵢ = sum((R Q⁻¹) ⊙ R)
  quad   <- sum((resids %*% mats$Q_inv) * resids)
  ll_base - 0.5 * quad
}

# ---- Numerical MLE for the coupled OU block ----
# Optimises (log a_p, log a_s, a_ps, a_sp, log σ_p, log σ_s).
# Input gains c_p = c_s ≡ 0: inference treats u(t) = 0 since the external
# input is never observed in practice. Dynamic variation is attributed entirely
# to the OU spectral structure.

ou_coupled_mle <- function(x_mat, a_p_init, a_s_init,
                           a_ps_init = 0, a_sp_init = 0, dt,
                           u_vec = NULL, c_p = 0, c_s = 0) {
  if (is.null(u_vec)) u_vec <- rep(0.0, nrow(x_mat))

  init_p <- ou_mle(x_mat[, 1L], u_vec, a_p_init, dt, c_fixed = c_p)
  init_s <- ou_mle(x_mat[, 2L], u_vec, a_s_init, dt, c_fixed = c_s)

  th0 <- c(log(init_p$a), log(init_s$a),
           max(a_ps_init, 1e-6), max(a_sp_init, 1e-6),
           log(max(init_p$sigma, 1e-5)), log(max(init_s$sigma, 1e-5)))

  neg_ll <- function(th) {
    a_p_v  <- exp(th[1L]);  a_s_v  <- exp(th[2L])
    a_ps_v <- th[3L];       a_sp_v <- th[4L]
    if (a_ps_v < 0 || a_sp_v < 0) return(1e10)
    ll <- tryCatch(
      ou_coupled_log_lik(a_ps_v, a_sp_v, a_p_v, a_s_v,
                         exp(th[5L]), exp(th[6L]),
                         c_p, c_s,
                         x_mat, u_vec, dt),
      error = function(e) -Inf
    )
    if (is.finite(ll)) -ll else 1e10
  }

  res <- optim(th0, neg_ll, method = "L-BFGS-B",
               lower = c(-Inf, -Inf, 0, 0, -Inf, -Inf),
               control = list(maxit = 500L, factr = 1e8))
  th   <- res$par

  H    <- numerical_hessian(neg_ll, th)
  V    <- tryCatch(solve(H), error = function(e) matrix(NA_real_, 6L, 6L))
  s_se <- function(k) sqrt(max(V[k, k], 0, na.rm = TRUE))

  a_p_h <- exp(th[1L]);  a_s_h <- exp(th[2L])
  sig_p <- exp(th[5L]);  sig_s <- exp(th[6L])
  list(
    a_p = a_p_h,    a_s = a_s_h,
    a_ps = th[3L],  a_sp = th[4L],
    sigma_p = sig_p, sigma_s = sig_s,
    c_p = c_p,       c_s = c_s,       # carry through the fixed known constants
    a_p_se     = a_p_h  * s_se(1L),
    a_s_se     = a_s_h  * s_se(2L),
    a_ps_se    = s_se(3L),
    a_sp_se    = s_se(4L),
    sigma_p_se = sig_p  * s_se(5L),
    sigma_s_se = sig_s  * s_se(6L),
    c_p_se = NA_real_,  c_s_se = NA_real_,
    log_lik     = if (res$value >= 1e9) NA_real_ else -res$value,
    convergence = res$convergence,
    n_obs       = nrow(x_mat) - 1L
  )
}

# ---- Block 3: Analytic IG observation MLE ----
#
# Given inter-beat intervals tau_k and known net-drive delta_k = s(t_k) - p(t_k),
# the conditional mean is mu_k = mu_0 * exp(-delta_k).
#
# Score equations yield closed-form solutions:
#   mu_0_hat = sum(tau_k * exp(2*delta_k)) / sum(exp(delta_k))
#   kappa_hat = N / sum((tau_k - mu_k_hat)^2 / (mu_k_hat^2 * tau_k))
#
# Derivation in supplementary. Reduces to standard IG MLE when delta_k = 0.

ig_obs_mle <- function(tau_vec, delta_vec) {
  stopifnot(length(tau_vec) == length(delta_vec), length(tau_vec) > 1L)

  if (any(tau_vec <= 0))
    stop("ig_obs_mle: tau_vec contains non-positive intervals — check spike times")
  if (any(!is.finite(delta_vec)))
    stop("ig_obs_mle: delta_vec contains non-finite values — check state trajectory")

  g       <- exp(-delta_vec)         # g_k = mu_k / mu_0  (=1 when delta=0)
  w       <- exp( delta_vec)         # w_k = exp(delta_k); direct avoids 1/g overflow

  # Closed-form MLE for mu_0
  denom_w <- sum(w)
  if (!is.finite(denom_w) || denom_w < .Machine$double.eps) {
    warning("ig_obs_mle: denominator sum(exp(delta)) near zero; returning NAs.")
    return(list(mu0=NA_real_, kappa=NA_real_, mu0_se=NA_real_, kappa_se=NA_real_,
                log_lik=NA_real_, n_beats=length(tau_vec)))
  }
  mu0_hat <- sum(tau_vec * w^2) / denom_w

  # Closed-form MLE for kappa (given mu_0_hat)
  mu_k    <- pmax(mu0_hat * g, .Machine$double.eps^0.5)
  psi     <- sum((tau_vec - mu_k)^2 / (mu_k^2 * tau_vec))
  n       <- length(tau_vec)
  kappa_hat <- n / max(psi, 1e-10)
  if (!is.finite(kappa_hat) || kappa_hat > 1e6) {
    warning(sprintf(
      "ig_obs_mle: kappa_hat = %.3g (capped at 1e6). psi = %.3g, n = %d.",
      kappa_hat, psi, n))
    kappa_hat <- 1e6
  }

  # Standard errors via analytic FIM (derived in supplementary):
  #   I(log mu_0) = kappa * sum(exp(delta_k)) / mu_0
  #   I(log kappa) = N / 2
  #   Cross term = 0 (block diagonal in log-space)
  i_log_mu0   <- kappa_hat * sum(w) / mu0_hat
  # Derivation: I(log mu_0) = kappa * sum(1/mu_k) = (kappa/mu_0) * sum(exp(delta_k))
  # where w_k = exp(delta_k) = mu_0/mu_k.  Reduces to N*kappa/mu_0 when delta=0 everywhere

  i_log_kappa <- n / 2
  # SEs are unreliable if kappa is at its cap; flag with NA
  mu0_se   <- if (kappa_hat >= 1e6) NA_real_ else mu0_hat   / sqrt(i_log_mu0)
  kappa_se <- if (kappa_hat >= 1e6) NA_real_ else kappa_hat / sqrt(i_log_kappa)

  ll <- sum(log_ig_pdf(tau_vec, mu_k, kappa_hat))

  rho_hat <- sqrt(mu0_hat / kappa_hat)
  # SE(rho): exact analytical formula via the (log mu_0, log rho) FIM.
  # At the MLE, d²L/d(log mu_0)d(log kappa) = 0 exactly (the score
  # d(log mu_0)L = -kappa*mu_0*dQ/dmu_0/2 vanishes, and differentiating
  # by log kappa leaves the same expression = 0). The Jacobian
  # d(log kappa)/d(log rho) = -2 then transforms the block-diagonal
  # (log mu_0, log kappa) FIM to:
  #   FIM(log mu_0, log rho) = [[i_log_mu0 + n/2,  -n],
  #                             [-n,              2*n]]
  # det = 2*n*i_log_mu0  (exact; no cancellation).
  # Var(log rho) = (i_log_mu0 + n/2) / (2*n*i_log_mu0).
  det_fim     <- 2 * n * i_log_mu0
  var_log_rho <- if (is.finite(det_fim) && det_fim > 0)
    (i_log_mu0 + n / 2) / det_fim else NA_real_
  rho_se      <- if (!is.na(var_log_rho)) rho_hat * sqrt(max(var_log_rho, 0)) else NA_real_


  list(mu0 = mu0_hat, kappa = max(kappa_hat, 1e-4),
       rho = rho_hat,
       mu0_se = mu0_se, kappa_se = kappa_se, rho_se = rho_se,
       log_lik = ll, n_beats = n)
}

# ---- Combined conditional MLE for all free parameters ----

full_conditional_mle <- function(sim_res, input_fn = NULL) {
  fp    <- sim_res$params$free
  dt    <- sim_res$time[2L] - sim_res$time[1L]

  # Resolve input: caller-supplied > stored in sim_res > zero.
  # c_p and c_s are treated as KNOWN fixed constants (not estimated).
  # Supplying the actual u(t) removes the input-driven drift from OU residuals,
  # sharpening identification of (a, σ) and eliminating the θ_true ≠ θ* mismatch.
  if (is.null(input_fn))
    input_fn <- if (!is.null(sim_res$input_fn)) sim_res$input_fn else function(t) 0
  u_vec <- vapply(sim_res$time, input_fn, numeric(1L))

  use_coupled <- (abs(fp$a_ps) + abs(fp$a_sp)) > 1e-12

  if (use_coupled) {
    x_mat  <- cbind(sim_res$p, sim_res$s)
    mle_ou <- ou_coupled_mle(x_mat,
                             a_p_init = fp$a_p, a_s_init = fp$a_s,
                             a_ps_init = fp$a_ps, a_sp_init = fp$a_sp, dt,
                             u_vec = u_vec, c_p = fp$c_p, c_s = fp$c_s)
    mle_p  <- list(a = mle_ou$a_p,   sigma = mle_ou$sigma_p, c_gain = mle_ou$c_p,
                   a_se = mle_ou$a_p_se, sigma_se = mle_ou$sigma_p_se,
                   c_se = mle_ou$c_p_se, log_lik = NA_real_)
    mle_s  <- list(a = mle_ou$a_s,   sigma = mle_ou$sigma_s, c_gain = mle_ou$c_s,
                   a_se = mle_ou$a_s_se, sigma_se = mle_ou$sigma_s_se,
                   c_se = mle_ou$c_s_se, log_lik = NA_real_)
  } else {
    mle_ou <- NULL
    mle_p  <- ou_mle(sim_res$p, u_vec, a_init = fp$a_p, dt, c_fixed = fp$c_p)
    mle_s  <- ou_mle(sim_res$s, u_vec, a_init = fp$a_s, dt, c_fixed = fp$c_s)
  }

  spk <- sim_res$spikes
  if (length(spk) < 2L) {
    na_out <- setNames(rep(NA_real_, 10),
                       c("a_p","a_s","sigma_p","sigma_s","mu0","rho","c_p","c_s","a_ps","a_sp"))
    return(c(as.list(na_out), list(n_beats = 0L)))
  }

  tau_vec   <- diff(spk)
  delta_v   <- compute_effective_delta(spk, sim_res$time, sim_res$delta)
  valid_idx <- is.finite(delta_v)

  if (sum(valid_idx) < 2L) {
    warning("full_conditional_mle: fewer than 2 valid IBIs after NA filtering.")

    return(c(as.list(setNames(rep(NA_real_, 14L),
                              c("a_p","a_s","a_ps","a_sp","sigma_p","sigma_s",
                                "c_p","c_s","mu0","kappa","rho","ll_obs","ll_total","n_beats"))),

             list(a_p_se = NA_real_, a_s_se = NA_real_, sigma_p_se = NA_real_,
                  sigma_s_se = NA_real_, mu0_se = NA_real_, rho_se = NA_real_,
                  a_ps_se = NA_real_, a_sp_se = NA_real_, c_p_se = NA_real_, c_s_se = NA_real_)))
  }
  mle_obs <- ig_obs_mle(tau_vec[valid_idx], delta_v[valid_idx])

  ll_ou <- if (!is.null(mle_ou)) mle_ou$log_lik
  else mle_p$log_lik + mle_s$log_lik

  list(
    a_p      = mle_p$a,      a_p_se      = mle_p$a_se,
    a_s      = mle_s$a,      a_s_se      = mle_s$a_se,
    a_ps     = if (!is.null(mle_ou)) mle_ou$a_ps    else fp$a_ps,
    a_sp     = if (!is.null(mle_ou)) mle_ou$a_sp    else fp$a_sp,
    a_ps_se  = if (!is.null(mle_ou)) mle_ou$a_ps_se else NA_real_,
    a_sp_se  = if (!is.null(mle_ou)) mle_ou$a_sp_se else NA_real_,
    sigma_p  = mle_p$sigma,  sigma_p_se  = mle_p$sigma_se,
    sigma_s  = mle_s$sigma,  sigma_s_se  = mle_s$sigma_se,
    c_p      = mle_p$c_gain, c_p_se      = mle_p$c_se,
    c_s      = mle_s$c_gain, c_s_se      = mle_s$c_se,
    mu0      = mle_obs$mu0,  mu0_se      = mle_obs$mu0_se,
    kappa    = mle_obs$kappa,
    rho      = mle_obs$rho,  rho_se      = mle_obs$rho_se,
    n_beats  = mle_obs$n_beats,
    ll_obs   = mle_obs$log_lik,
    ll_total = ll_ou + mle_obs$log_lik
  )
}

# ---- Moment consistency check ----
# Compares observed HRV statistics against the values predicted by the
# SDE-IG model at the MLE. Provides a pre-submission specification test
# complementary to the time-rescaling diagnostic.
#
# Under correct model specification, each relative residual should be
# small (< 5% in absolute value). Systematic departures indicate either
# a local-optimum convergence failure or genuine model misspecification.
#
# Arguments
#   rr_vec  : raw RR interval vector (seconds)
#   mle     : output of full_conditional_mle() or a list with named
#             fields a_p, a_s, sigma_p, sigma_s, mu0, rho
#
# Returns a data.frame suitable for knitr::kable().

check_moment_consistency <- function(rr_vec, mle) {
  stopifnot(is.numeric(rr_vec), length(rr_vec) > 10L, all(rr_vec > 0))

  # Normalise field names: full_conditional_mle() uses 'mu0';
  # pp_mle free_hat uses 'mu_0'.  Accept both.
  mu0_val <- if (!is.null(mle$mu0))  mle$mu0  else
    if (!is.null(mle$mu_0)) mle$mu_0 else
      stop("check_moment_consistency: mle must contain 'mu0' or 'mu_0'")

  pred <- predicted_hrv_moments(
    a_p     = mle$a_p,
    a_s     = mle$a_s,
    sigma_p = mle$sigma_p,
    sigma_s = mle$sigma_s,
    mu_0    = mu0_val,
    rho     = mle$rho
  )

  # ---- Observed statistics ----
  n       <- length(rr_vec)
  e_obs   <- mean(rr_vec)
  sd_obs  <- sd(rr_vec)                          # SDNN
  rmssd_obs <- sqrt(mean(diff(rr_vec)^2))

  # Spectral band powers via Lomb-Scargle on the tachogram.
  # Beat times are reconstructed assuming the first beat at t = 0.
  beat_t  <- cumsum(c(0, rr_vec))[-1L]           # midpoints would be more correct
  # Fast approximation: FFT on the evenly resampled tachogram (4 Hz)
  dt_rs   <- 0.25
  t_rs    <- seq(beat_t[1L], tail(beat_t, 1L), by = dt_rs)
  rr_rs   <- approx(beat_t, rr_vec, xout = t_rs, rule = 2)$y
  n_rs    <- length(rr_rs)
  n_half  <- floor(n_rs / 2L) + 1L
  freq_rs <- (seq_len(n_half) - 1L) / (n_rs * dt_rs)   # exact, uniform spacing
  df_rs   <- 1.0 / (n_rs * dt_rs)                       # exact frequency resolution

  psd_rs  <- (Mod(fft(rr_rs - mean(rr_rs)))[seq_len(n_half)])^2 /
    (n_rs * (1.0 / dt_rs))

  # Convert to one-sided PSD.
  # Even n_rs: last bin is Nyquist — double all bins 2:(n_half-1).
  # Odd n_rs:  no Nyquist bin — double all bins 2:n_half.
  if (n_half > 2L) {
    last_to_double <- if (n_rs %% 2L == 0L) n_half - 1L else n_half
    psd_rs[2L:last_to_double] <- 2.0 * psd_rs[2L:last_to_double]
  }
  idx_lf  <- freq_rs >= 0.04 & freq_rs <= 0.15
  idx_hf  <- freq_rs >= 0.15 & freq_rs <= 0.40
  lf_obs  <- if (any(idx_lf)) sum(psd_rs[idx_lf]) * df_rs else NA_real_
  hf_obs  <- if (any(idx_hf)) sum(psd_rs[idx_hf]) * df_rs else NA_real_
  lfhf_obs <- lf_obs / max(hf_obs, 1e-12)

  # ---- Residuals ----
  rel_err <- function(obs, pred_v)
    ifelse(abs(obs) > 1e-10, 100 * (obs - pred_v) / abs(obs), NA_real_)

  data.frame(
    Statistic = c(
      "Mean RR (s)",
      "SDNN (s)",
      "RMSSD (s)",
      "LF power (s\\textsuperscript{2})",
      "HF power (s\\textsuperscript{2})",
      "LF/HF ratio"
    ),
    Observed  = sprintf("%.4f", c(e_obs, sd_obs, rmssd_obs,
                                  lf_obs, hf_obs, lfhf_obs)),
    Predicted = sprintf("%.4f", c(pred$e_tau, pred$sdnn, pred$rmssd,
                                  pred$lf_power, pred$hf_power, pred$lf_hf)),
    Rel_error = sprintf("%.1f\\%%", c(
      rel_err(e_obs,    pred$e_tau),
      rel_err(sd_obs,   pred$sdnn),
      rel_err(rmssd_obs, pred$rmssd),
      rel_err(lf_obs,   pred$lf_power),
      rel_err(hf_obs,   pred$hf_power),
      rel_err(lfhf_obs, pred$lf_hf)
    )),
    stringsAsFactors = FALSE
  )
}

# ---- Observed information matrix (negative Hessian) for IG block ----
#
# Parameterization: theta = (log mu_0, log kappa) for numerical stability.
# Returns the 2x2 FIM, its inverse (parameter covariance), SE, and
# condition number.

ig_obs_fim <- function(mu0, rho_val, tau_vec, delta_vec) {
  # Analytical observed FIM in (log mu_0, log rho) parameterisation.
  # The (log mu_0, log kappa) observed FIM is block-diagonal at the MLE
  # (cross-term exactly zero; see ig_obs_mle comment for proof).
  # Jacobian d(log kappa)/d(log rho) = -2 gives:
  #   FIM = [[i_mu0 + n/2,  -n],
  #          [-n,         2*n]]
  # with i_mu0 = kappa * sum(exp(delta_k)) / mu_0.
  # det = 2*n*i_mu0  (exact).
  # FIM^{-1} = [[2n, n], [n, i_mu0 + n/2]] / det.
  n     <- length(tau_vec)
  kap   <- kappa_from_rho(mu0, rho_val)
  w     <- exp(delta_vec)
  i_mu0 <- kap * sum(w) / mu0
  fim   <- matrix(c(i_mu0 + n / 2, -n, -n, 2 * n), 2L, 2L)
  det_v <- 2 * n * i_mu0
  fim_inv <- if (is.finite(det_v) && det_v > 0)
    matrix(c(2 * n, n, n, i_mu0 + n / 2), 2L, 2L) / det_v
  else
    matrix(NA_real_, 2L, 2L)
  eig_v <- tryCatch(
    eigen(fim, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, 2L)
  )
  list(
    fim  = fim,
    se   = sqrt(pmax(diag(fim_inv), 0, na.rm = TRUE)),
    cond = if (all(is.finite(eig_v)) && min(eig_v) > 0)
      max(eig_v) / min(eig_v) else NA_real_,
    eig  = eig_v
  )
}

# ---- Observed information matrix (negative Hessian) for one OU branch ----
#
# Parameterization: theta = (log a, log sigma) — 2D.
# c_gain is excluded: inference always uses u(t) ≡ 0, so c is unidentified
# and must not appear in the FIM (it would produce a structural zero eigenvalue).

ou_fim <- function(a_val, sigma_val, x_vec, dt, u_vec = NULL, c_fixed = 0) {
  th0 <- c(log(a_val), log(sigma_val))
  n   <- length(x_vec) - 1L
  idx <- seq_len(n)
  if (is.null(u_vec)) u_vec <- rep(0.0, length(x_vec))

  nll_2d <- function(th) {
    a_v   <- exp(th[1L]); sig_v <- exp(th[2L])
    if (a_v <= 0 || sig_v <= 0) return(1e10)
    e_adt <- exp(-a_v * dt)
    z_v   <- u_vec[idx] / a_v * (-expm1(-a_v * dt))
    y_v   <- x_vec[-1L] - x_vec[idx] * e_adt - c_fixed * z_v  # remove known input
    C_v   <- (-expm1(-2 * a_v * dt)) / (2 * a_v)
    ll_v  <- tryCatch(sum(dnorm(y_v, 0, sig_v * sqrt(C_v), log = TRUE)),
                      error = function(e) -Inf)
    if (is.finite(ll_v)) -ll_v else 1e10
  }

  H     <- numerical_hessian(nll_2d, th0)
  det_H <- H[1L,1L] * H[2L,2L] - H[1L,2L] * H[2L,1L]
  fim_inv <- if (is.finite(det_H) && det_H > 0) {
    matrix(c(H[2L,2L], -H[2L,1L], -H[1L,2L], H[1L,1L]), 2L, 2L) / det_H
  } else {
    matrix(NA_real_, 2L, 2L)
  }

  list(
    fim  = H,
    se   = sqrt(abs(diag(fim_inv))),   # SE(log a), SE(log sigma)
    cond = tryCatch(kappa(H), error = function(e) NA_real_),
    eig  = tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                    error = function(e) rep(NA_real_, 2L))
  )
}

# ---- Full 6x6 block-diagonal FIM ----
#
# Parameterization: (log sigma_p, c_p, log sigma_s, c_s, log mu_0, log kappa)
# Returns the full FIM, all eigenvalues, overall condition number,
# and per-parameter SEs.

# Internal: full numerical FIM for the coupled case.
# Parameterisation: (b_ps, b_sp, log σ_p, log σ_s, c_p, c_s, log μ₀, log κ).
# The block-diagonal structure no longer holds when b_ps or b_sp ≠ 0.

.full_fim_coupled <- function(sim_res, mle, dt, u_vec = NULL) {
  x_mat   <- cbind(sim_res$p, sim_res$s)
  if (is.null(u_vec)) u_vec <- rep(0.0, nrow(x_mat))
  fp_sim  <- sim_res$params$free
  spk     <- sim_res$spikes
  tau_vec <- diff(spk)
  delta_v <- compute_effective_delta(spk, sim_res$time, sim_res$delta)

  # Mirror the NA-filtering of full_conditional_fim (uncoupled path):
  # IBIs shorter than dt produce NA from compute_effective_delta; discard them.
  valid_mask <- is.finite(delta_v)
  tau_vec    <- tau_vec[valid_mask]
  delta_v    <- delta_v[valid_mask]

  # 8 free parameters: (log a_p, log a_s, a_ps, a_sp,
  #                      log σ_p, log σ_s, log μ₀, log ρ)
  # c_p = c_s = 0 throughout (not estimated from RR data).
  th0 <- c(log(mle$a_p),  log(mle$a_s),
           mle$a_ps,       mle$a_sp,
           log(mle$sigma_p), log(mle$sigma_s),
           log(mle$mu0),  log(mle$rho))

  obj <- function(th) {
    a_p_v <- exp(th[1L]);  a_s_v  <- exp(th[2L])
    a_ps_v <- th[3L];      a_sp_v <- th[4L]
    if (a_ps_v < 0 || a_sp_v < 0 || a_p_v * a_s_v - a_ps_v * a_sp_v <= 0)
      return(-Inf)

    ll_ou <- tryCatch(
      ou_coupled_log_lik(a_ps_v, a_sp_v, a_p_v, a_s_v,
                         exp(th[5L]), exp(th[6L]),
                         fp_sim$c_p, fp_sim$c_s,
                         x_mat, u_vec, dt),
      error = function(e) -Inf)
    mu_k   <- exp(th[7L]) * exp(-delta_v)
    kap_v  <- kappa_from_rho(exp(th[7L]), exp(th[8L]))
    ll_obs <- sum(log_ig_pdf(tau_vec, mu_k, kap_v))
    ll_ou + ll_obs
  }

  # obj returns log-likelihood; Hessian of LL = -FIM, so negate.
  FIM <- -numerical_hessian(obj, th0)

  pnames <- c("log(a_p)", "log(a_s)", "a_ps", "a_sp",
              "log(sigma_p)", "log(sigma_s)", "log(mu_0)", "log(rho)")
  dimnames(FIM) <- list(pnames, pnames)
  eig_all <- tryCatch(
    eigen(FIM, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(FIM)))
  list(FIM         = FIM,
       eigenvalues = eig_all,
       cond_number = max(eig_all) / min(eig_all[eig_all > 0]),
       param_names = pnames,
       blocks      = list(p   = list(eig = NA, cond = NA),
                          s   = list(eig = NA, cond = NA),
                          obs = list(eig = tail(eig_all, 2L), cond = NA)))
}

# ---- Full 6x6 block-diagonal FIM ----
#
# Parameterization: (log a_p, log sigma_p, log a_s, log sigma_s, log mu_0, log rho)
# c_p and c_s are excluded: they are simulation-only and unidentified from u(t) ≡ 0.
# Returns the full FIM, all eigenvalues, overall condition number, and per-param SEs.

full_conditional_fim <- function(sim_res, input_fn = NULL) {
  dt  <- sim_res$time[2L] - sim_res$time[1L]
  if (is.null(input_fn)) {
    if (!is.null(sim_res$input_fn)) {
      input_fn <- sim_res$input_fn
    } else {
      input_fn <- function(t) { 0 }
    }
  }
  u_vec  <- vapply(sim_res$time, input_fn, numeric(1L))
  fp_sim <- sim_res$params$free
  mle    <- full_conditional_mle(sim_res, input_fn = input_fn)

  if ((abs(fp_sim$a_ps) + abs(fp_sim$a_sp)) > 1e-12) {
    return(.full_fim_coupled(sim_res, mle, dt, u_vec = u_vec))
  }

  fim_p <- ou_fim(mle$a_p, mle$sigma_p, sim_res$p, dt,
                  u_vec = u_vec, c_fixed = fp_sim$c_p)
  fim_s <- ou_fim(mle$a_s, mle$sigma_s, sim_res$s, dt,
                  u_vec = u_vec, c_fixed = fp_sim$c_s)

  spk        <- sim_res$spikes
  tau_vec    <- diff(spk)
  delta_v    <- compute_effective_delta(spk, sim_res$time, sim_res$delta)
  valid_mask <- is.finite(delta_v)
  fim_obs    <- ig_obs_fim(mle$mu0, mle$rho,
                           tau_vec[valid_mask], delta_v[valid_mask])

  # 6×6 block-diagonal FIM in:
  # (log a_p, log sigma_p, log a_s, log sigma_s, log mu_0, log rho)
  FIM <- matrix(0, 6, 6)
  FIM[1:2, 1:2] <- fim_p$fim
  FIM[3:4, 3:4] <- fim_s$fim
  FIM[5:6, 5:6] <- fim_obs$fim

  param_names <- c("log(a_p)", "log(sigma_p)",
                   "log(a_s)", "log(sigma_s)",
                   "log(mu_0)", "log(rho)")
  dimnames(FIM) <- list(param_names, param_names)

  eig_all <- tryCatch(
    eigen(FIM, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, 6)
  )

  list(
    FIM         = FIM,
    eigenvalues = eig_all,
    cond_number = max(eig_all) / min(eig_all[eig_all > 0]),
    param_names = param_names,
    blocks      = list(p   = fim_p,
                       s   = fim_s,
                       obs = fim_obs)
  )
}

# ---- Profile likelihood for each free parameter ----
#
# For each parameter, the profile is computed by maximizing over the
# remaining parameters in its block.  Exception: profiles of sigma_p,
# sigma_s, c_p, c_s hold a_p / a_s at the simulation reference value
# (fp$a_p, fp$a_s) rather than profiling over them jointly.  This is a
# close approximation because ou_mle() recovers the decay rates with <1%
# relative bias, making fp$a_p ≈ mle$a_p.

profile_lik_one <- function(param, grid, sim_res, mle, input_fn = NULL) {
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  # c_p and c_s are fixed known constants.  u(t) is supplied so their
  # contribution is removed from OU residuals before profiling (a, σ).
  if (is.null(input_fn))
    input_fn <- if (!is.null(sim_res$input_fn)) sim_res$input_fn else function(t) 0
  u_vec <- vapply(sim_res$time, input_fn, numeric(1L))
  fp    <- sim_res$params$free

  a_p_mle <- mle$a_p
  a_s_mle <- mle$a_s

  # Branch pre-computations: residuals after removing known input contribution.
  # RSS and C_a are constant w.r.t. sigma, enabling an O(1) analytical formula
  # per sigma-profile grid point instead of an O(n_time) ou_log_lik_at call.
  n_p    <- length(sim_res$p) - 1L
  idx_p  <- seq_len(n_p)
  e_p    <- exp(-a_p_mle * dt)
  z_p    <- u_vec[idx_p] / a_p_mle * (-expm1(-a_p_mle * dt))
  y_p    <- sim_res$p[-1L] - sim_res$p[idx_p] * e_p - fp$c_p * z_p
  C_ap   <- (-expm1(-2 * a_p_mle * dt)) / (2 * a_p_mle)
  rss_p  <- sum(y_p^2)
  ll_const_p <- -n_p / 2 * (log(2 * pi) + log(C_ap))

  n_s    <- length(sim_res$s) - 1L
  idx_s  <- seq_len(n_s)
  e_s    <- exp(-a_s_mle * dt)
  z_s    <- u_vec[idx_s] / a_s_mle * (-expm1(-a_s_mle * dt))
  y_s    <- sim_res$s[-1L] - sim_res$s[idx_s] * e_s - fp$c_s * z_s
  C_as   <- (-expm1(-2 * a_s_mle * dt)) / (2 * a_s_mle)
  rss_s  <- sum(y_s^2)
  ll_const_s <- -n_s / 2 * (log(2 * pi) + log(C_as))

  tau_vec <- diff(sim_res$spikes)
  delta_v <- compute_effective_delta(sim_res$spikes, sim_res$time, sim_res$delta)

  # Mirror the NA-filtering of full_conditional_mle: discard IBIs whose
  # interval-averaged drive could not be computed (IBI < dt).
  valid_mask <- is.finite(delta_v)
  tau_vec    <- tau_vec[valid_mask]
  delta_v    <- delta_v[valid_mask]
  if (length(tau_vec) < 2L) {
    return(rep(-Inf, length(grid)))
  }

  g       <- exp(-delta_v)
  w       <- exp(delta_v)
  mu0_hat <- sum(tau_vec * w^2) / sum(w)

  # Pre-compute rho-profile constants (O(n_beats), used for every grid point).
  rho_C1   <- sum(tau_vec * w^2)    # = sum(tau_k * exp(2*delta_k))
  rho_C2   <- sum(1 / tau_vec)
  rho_nobs <- length(tau_vec)

  # Hoist coupled-model matrix: O(n_time) allocation, constant across all grid points.
  x_mat_coupled <- if ((abs(sim_res$params$free$a_ps) +
                        abs(sim_res$params$free$a_sp)) > 1e-12)
    cbind(sim_res$p, sim_res$s) else NULL

  # profile_at_a_internal: profile over a with c_known fixed and σ eliminated.
  profile_at_a_internal <- function(x_vec_b, u_vec_b, c_known, a_v, dt_b) {
    n_b   <- length(x_vec_b) - 1L;  idx_b <- seq_len(n_b)
    e_adt <- exp(-a_v * dt_b)
    z_b   <- u_vec_b[idx_b] / a_v * (-expm1(-a_v * dt_b))
    y_b   <- x_vec_b[-1L] - x_vec_b[idx_b] * e_adt - c_known * z_b
    C_b   <- (-expm1(-2 * a_v * dt_b)) / (2 * a_v)
    sig2  <- max(mean(y_b^2) / C_b, 1e-14)
    list(ll = sum(dnorm(y_b, 0, sqrt(sig2 * C_b), log = TRUE)))
  }

  # ---- Coupled parameters: warm-start sequential loop ----
  # Pre-compute x_mat ONCE outside the loop (not inside vapply per grid point).
  # Warm starts: pass previous converged solution as initialization for next point.
  if (param %in% c("a_ps", "a_sp")) {
    x_mat  <- cbind(sim_res$p, sim_res$s)
    is_ps  <- (param == "a_ps")

    # Initialize from the joint MLE; th_warm is updated after each convergence.
    # Initialize from the joint MLE; th_warm is updated after each convergence.
    th_warm <- if (is_ps) {
      c(log(mle$a_p), log(mle$a_s), max(mle$a_sp, 1e-6),
        log(max(mle$sigma_p, 1e-5)), log(max(mle$sigma_s, 1e-5)))
    } else {
      c(log(mle$a_p), log(mle$a_s), max(mle$a_ps, 1e-6),
        log(max(mle$sigma_p, 1e-5)), log(max(mle$sigma_s, 1e-5)))
    }
    # Fixed cold-start reference: joint MLE, never mutated.
    # Used alongside the warm start to prevent cascading of local optima
    # that produce the sawtooth pattern seen in jagged profile curves.
    th_cold <- th_warm

    ll_out <- numeric(length(grid))
    for (k in seq_along(grid)) {
      v <- grid[k]
      if (!is.finite(v) || v < 0) { ll_out[k] <- -Inf; next }

      nll_k <- if (is_ps) {
        local({
          vv <- v
          function(th) {
            a_p_v <- exp(th[1L]);  a_s_v <- exp(th[2L]);  a_sp_v <- th[3L]
            if (a_sp_v < 0 || a_p_v * a_s_v - vv * a_sp_v <= 0) return(1e10)
            ll <- tryCatch(
              ou_coupled_log_lik(vv, a_sp_v, a_p_v, a_s_v,
                                 exp(th[4L]), exp(th[5L]),
                                 fp$c_p, fp$c_s,
                                 x_mat, u_vec, dt),
              error = function(e) -Inf)
            if (is.finite(ll)) -ll else 1e10
          }
        })
      } else {
        local({
          vv <- v
          function(th) {
            a_p_v <- exp(th[1L]);  a_s_v <- exp(th[2L]);  a_ps_v <- th[3L]
            if (a_ps_v < 0 || a_p_v * a_s_v - a_ps_v * vv <= 0) return(1e10)
            ll <- tryCatch(
              ou_coupled_log_lik(a_ps_v, vv, a_p_v, a_s_v,
                                 exp(th[4L]), exp(th[5L]),
                                 fp$c_p, fp$c_s,
                                 x_mat, u_vec, dt),
              error = function(e) -Inf)
            if (is.finite(ll)) -ll else 1e10
          }
        })
      }

      # Warm-start optimization: follows the profile curve efficiently when smooth.
      res_w <- tryCatch(
        optim(th_warm, nll_k, method = "L-BFGS-B",
              lower = c(-Inf, -Inf, 0, -Inf, -Inf),
              control = list(maxit = 300L, factr = 1e7)),
        error = function(e) list(value = 1e10, par = th_warm, convergence = 99L)
      )
      # Cold-start from joint MLE: independent fallback that cannot inherit
      # a local-minimum cascade from a previous grid point.
      res_c <- tryCatch(
        optim(th_cold, nll_k, method = "L-BFGS-B",
              lower = c(-Inf, -Inf, 0, -Inf, -Inf),
              control = list(maxit = 300L, factr = 1e7)),
        error = function(e) list(value = 1e10, par = th_cold, convergence = 99L)
      )
      res <- if (res_c$value < res_w$value) res_c else res_w

      ll_out[k] <- if (is.finite(-res$value)) -res$value else -Inf
      # Only advance the warm start from non-degenerate solutions.
      if (is.finite(res$value) && res$value < 1e9) th_warm <- res$par
    }
    return(ll_out)
  }

  # ---- All other parameters: standard vapply ----
  vapply(grid, function(v) {
    if (!is.finite(v)) return(-Inf)
    switch(param,
           mu0 = {
             if (v <= 0) return(-Inf)
             mu_k <- v * g
             psi  <- sum((tau_vec - mu_k)^2 / (mu_k^2 * tau_vec))
             k    <- min(length(tau_vec) / max(psi, 1e-10), 1e6)
             sum(log_ig_pdf(tau_vec, mu_k, k))
           },
           sigma_p = {
             if (v <= 0) return(-Inf)
             use_coupled_loc <- (abs(sim_res$params$free$a_ps) +
                                   abs(sim_res$params$free$a_sp)) > 1e-12
             if (!use_coupled_loc) {
               # Exact analytical formula: constant in a; O(1) per grid point.
               ll_const_p - n_p * log(v) - rss_p / (2 * v^2 * C_ap)
             } else {
               x_mat_loc <- x_mat_coupled
               th0_loc <- c(log(mle$a_p), log(mle$a_s),
                            max(mle$a_ps, 1e-6), max(mle$a_sp, 1e-6),
                            log(max(mle$sigma_s, 1e-5)))
               local({
                 vv <- v
                 nll_loc <- function(th) {
                   a_p_v  <- exp(th[1L]); a_s_v  <- exp(th[2L])
                   a_ps_v <- th[3L];      a_sp_v <- th[4L]
                   sig_s_v <- exp(th[5L])
                   if (a_ps_v < 0 || a_sp_v < 0 ||
                       a_p_v * a_s_v - a_ps_v * a_sp_v <= 0) return(1e10)
                   ll <- tryCatch(
                     ou_coupled_log_lik(a_ps_v, a_sp_v, a_p_v, a_s_v,
                                        vv, sig_s_v,
                                        fp$c_p, fp$c_s, x_mat_loc, u_vec, dt),
                     error = function(e) -Inf)
                   if (is.finite(ll)) -ll else 1e10
                 }
                 res_loc <- tryCatch(
                   optim(th0_loc, nll_loc, method = "L-BFGS-B",
                         lower = c(-Inf, -Inf, 0, 0, -Inf),
                         control = list(maxit = 200L, factr = 1e8)),
                   error = function(e) list(value = 1e10))
                 if (is.finite(-res_loc$value)) -res_loc$value else -Inf
               })
             }
           },
           sigma_s = {
             if (v <= 0) return(-Inf)
             use_coupled_loc <- (abs(sim_res$params$free$a_ps) +
                                   abs(sim_res$params$free$a_sp)) > 1e-12
             if (!use_coupled_loc) {
               # Exact analytical formula: constant in a; O(1) per grid point.
               ll_const_s - n_s * log(v) - rss_s / (2 * v^2 * C_as)
             } else {
               x_mat_loc <- x_mat_coupled
               th0_loc <- c(log(mle$a_p), log(mle$a_s),
                            max(mle$a_ps, 1e-6), max(mle$a_sp, 1e-6),
                            log(max(mle$sigma_p, 1e-5)))
               local({
                 vv <- v
                 nll_loc <- function(th) {
                   a_p_v  <- exp(th[1L]); a_s_v  <- exp(th[2L])
                   a_ps_v <- th[3L];      a_sp_v <- th[4L]
                   sig_p_v <- exp(th[5L])
                   if (a_ps_v < 0 || a_sp_v < 0 ||
                       a_p_v * a_s_v - a_ps_v * a_sp_v <= 0) return(1e10)
                   ll <- tryCatch(
                     ou_coupled_log_lik(a_ps_v, a_sp_v, a_p_v, a_s_v,
                                        sig_p_v, vv,
                                        fp$c_p, fp$c_s, x_mat_loc, u_vec, dt),
                     error = function(e) -Inf)
                   if (is.finite(ll)) -ll else 1e10
                 }
                 res_loc <- tryCatch(
                   optim(th0_loc, nll_loc, method = "L-BFGS-B",
                         lower = c(-Inf, -Inf, 0, 0, -Inf),
                         control = list(maxit = 200L, factr = 1e8)),
                   error = function(e) list(value = 1e10))
                 if (is.finite(-res_loc$value)) -res_loc$value else -Inf
               })
             }
           },
           a_p = {
             if (v <= 0) return(-Inf)
             use_coupled_loc <- (abs(sim_res$params$free$a_ps) +
                                   abs(sim_res$params$free$a_sp)) > 1e-12
             if (!use_coupled_loc) {
               profile_at_a_internal(sim_res$p, u_vec, fp$c_p, v, dt)$ll
             } else {
               # In coupled model, profile a_p within the full joint OU LL,
               # maximising over (a_s, a_ps, a_sp, sigma_p, sigma_s).
               x_mat_loc  <- x_mat_coupled
               th0_loc <- c(log(mle$a_s),
                            max(mle$a_ps, 1e-6), max(mle$a_sp, 1e-6),
                            log(max(mle$sigma_p, 1e-5)),
                            log(max(mle$sigma_s, 1e-5)))
               local({
                 vv <- v
                 nll_loc <- function(th) {
                   a_s_v  <- exp(th[1L])
                   a_ps_v <- th[2L];  a_sp_v <- th[3L]
                   if (a_ps_v < 0 || a_sp_v < 0 ||
                       vv * a_s_v - a_ps_v * a_sp_v <= 0) return(1e10)
                   ll <- tryCatch(
                     ou_coupled_log_lik(a_ps_v, a_sp_v, vv, a_s_v,
                                        exp(th[4L]), exp(th[5L]),
                                        fp$c_p, fp$c_s,
                                        x_mat_loc, u_vec, dt),
                     error = function(e) -Inf)
                   if (is.finite(ll)) -ll else 1e10
                 }
                 -optim(th0_loc, nll_loc, method = "L-BFGS-B",
                        lower = c(-Inf, 0, 0, -Inf, -Inf),
                        control = list(maxit = 200L, factr = 1e8))$value
               })
             }
           },
           a_s = {
             if (v <= 0) return(-Inf)
             use_coupled_loc <- (abs(sim_res$params$free$a_ps) +
                                   abs(sim_res$params$free$a_sp)) > 1e-12
             if (!use_coupled_loc) {
               profile_at_a_internal(sim_res$s, u_vec, fp$c_s, v, dt)$ll
             } else {
               x_mat_loc  <- x_mat_coupled
               th0_loc <- c(log(mle$a_p),
                            max(mle$a_ps, 1e-6), max(mle$a_sp, 1e-6),
                            log(max(mle$sigma_p, 1e-5)),
                            log(max(mle$sigma_s, 1e-5)))
               local({
                 vv <- v
                 nll_loc <- function(th) {
                   a_p_v  <- exp(th[1L])
                   a_ps_v <- th[2L];  a_sp_v <- th[3L]
                   if (a_ps_v < 0 || a_sp_v < 0 ||
                       a_p_v * vv - a_ps_v * a_sp_v <= 0) return(1e10)
                   ll <- tryCatch(
                     ou_coupled_log_lik(a_ps_v, a_sp_v, a_p_v, vv,
                                        exp(th[4L]), exp(th[5L]),
                                        fp$c_p, fp$c_s,
                                        x_mat_loc, u_vec, dt),
                     error = function(e) -Inf)
                   if (is.finite(ll)) -ll else 1e10
                 }
                 -optim(th0_loc, nll_loc, method = "L-BFGS-B",
                        lower = c(-Inf, 0, 0, -Inf, -Inf),
                        control = list(maxit = 200L, factr = 1e8))$value
               })
             }
           },
           rho = {
             if (v <= 0) return(-Inf)
             disc    <- sqrt((rho_nobs * v^2)^2 + 4 * rho_C2 * rho_C1)
             mu0_opt <- (rho_nobs * v^2 + disc) / (2 * rho_C2)
             if (!is.finite(mu0_opt) || mu0_opt <= 0) return(-Inf)
             kap_opt <- kappa_from_rho(mu0_opt, v)
             sum(log_ig_pdf(tau_vec, mu0_opt * g, kap_opt))
           },
           -Inf
    )
  }, numeric(1L))
}

# ---- All profile likelihoods in one call ----

all_profile_likelihoods <- function(sim_res, n_grid = 60, width = 3.0,
                                    input_fn = NULL) {
  mle <- full_conditional_mle(sim_res, input_fn = input_fn)

  fp_sim      <- sim_res$params$free
  use_coupled <- (abs(fp_sim$a_ps) + abs(fp_sim$a_sp)) > 1e-12

  specs <- list(
    a_p     = list(center = mle$a_p,     log_scale = TRUE,
                   label = expression(a[p] ~ "(Hz)")),
    a_s     = list(center = mle$a_s,     log_scale = TRUE,
                   label = expression(a[s] ~ "(Hz)")),
    sigma_p = list(center = mle$sigma_p, log_scale = TRUE,
                   label = expression(sigma[p])),
    sigma_s = list(center = mle$sigma_s, log_scale = TRUE,
                   label = expression(sigma[s])),
    mu0     = list(center = mle$mu0,     log_scale = TRUE,
                   label = expression(mu[0] ~ "(s)")),
    rho     = list(center = mle$rho,     log_scale = TRUE,
                   label = expression(rho ~ "(baseline CV)"))
  )

  if (use_coupled) {
    specs$a_ps <- list(center = max(mle$a_ps, 1e-4), log_scale = TRUE,
                       width_override = 3.0, label = expression(a[ps]))
    specs$a_sp <- list(center = max(mle$a_sp, 1e-4), log_scale = TRUE,
                       width_override = 3.0, label = expression(a[sp]))
  }

  lapply(names(specs), function(p) {
    s <- specs[[p]]
    eff_width <- if (!is.null(s$width_override)) s$width_override else width
    grid <- if (s$log_scale) {
      exp(seq(log(s$center / eff_width), log(s$center * eff_width), length.out = n_grid))
    } else {
      spread <- max(abs(s$center) * eff_width, 0.5)
      seq(s$center - spread, s$center + spread, length.out = n_grid)
    }
    ll <- profile_lik_one(p, grid, sim_res, mle = mle, input_fn = input_fn)
    list(param = p, grid = grid, ll = ll,
         mle = s$center, label = s$label,
         ll_mle = max(ll, na.rm = TRUE))
  }) |> setNames(names(specs))
}
