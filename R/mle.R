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
  vapply(seq_along(spk[-1L]), function(k) {
    t_start <- spk[k]
    t_end   <- spk[k + 1L]
    idx <- which(time_vec > t_start & time_vec <= t_end)
    if (length(idx) > 0L) {
      -log(mean(exp(-delta_vec[idx])))
    } else {
      # No grid points inside (t_start, t_end] — IBI < dt.
      # This should not occur at dt = 0.005 s with physiological IBIs.
      # Using nearest-endpoint delta is equivalent to the old biased estimate;
      # flag it so the caller can investigate.
      warning(sprintf(
        "compute_effective_delta: no grid points in interval [%.4f, %.4f] (IBI=%.4f < dt). Using endpoint delta.",
        t_start, t_end, t_end - t_start))
      delta_vec[which.min(abs(time_vec - t_end))]
    }
  }, numeric(1L))
}

# Helper: reconstruct input vector from a sim_res object.
# Called once per analysis function; result cached by caller.
extract_u_vec <- function(sim_res) {
  if (!is.null(sim_res$input_fn))
    vapply(sim_res$time, sim_res$input_fn, numeric(1L))
  else
    rep(0.0, length(sim_res$time))
}

# ---- Block 1 / 2: Analytic OU MLE ----
#
# Model: x[i+1] | x[i] ~ N(x[i]*exp(-a*dt) + (c*u[i]/a)*(1-exp(-a*dt)),
#                           sigma^2/(2a)*(1-exp(-2a*dt)))
#
# Closed-form solution from OLS (for c) + method-of-moments (for sigma).

ou_mle <- function(x_vec, u_vec, a, dt) {
  n     <- length(x_vec) - 1L
  if (n < 2L) return(list(sigma = NA_real_, c_gain = NA_real_,
                           sigma_se = NA_real_, c_se = NA_real_,
                           log_lik = NA_real_))

  e_adt <- exp(-a * dt)
  z     <- u_vec[seq_len(n)] / a * (-expm1(-a * dt))  # input basis vector
  y     <- x_vec[-1L] - x_vec[seq_len(n)] * e_adt     # mean-subtracted increments

  # OLS for c_gain (exact MLE for linear Gaussian model)
  ssz   <- sum(z^2)
  c_hat <- if (ssz > 1e-14) sum(y * z) / ssz else 0
  resid <- y - c_hat * z

  # MLE for sigma (exact closed form)
  C_a      <- (-expm1(-2 * a * dt)) / (2 * a)   # Var(resid) / sigma^2
  sigma2   <- max(mean(resid^2) / C_a, 1e-14)
  sigma_hat <- sqrt(sigma2)

  # Standard errors (from observed information at MLE)
  #   SE(log sigma) = 1/sqrt(2*n)  =>  SE(sigma) = sigma/sqrt(2*n)
  #   SE(c) = sigma * sqrt(C_a / ssz)  [from OLS, homoscedastic]
  sigma_se <- sigma_hat / sqrt(2 * n)
  c_se     <- if (ssz > 1e-14) sigma_hat * sqrt(C_a / ssz) else Inf

  # Log-likelihood at MLE
  ll <- sum(dnorm(resid, mean = 0, sd = sigma_hat * sqrt(C_a), log = TRUE))

  list(sigma = sigma_hat, c_gain = c_hat,
       sigma_se = sigma_se, c_se = c_se,
       log_lik = ll, n_obs = n)
}

# ---- Bivariate coupled-OU log-likelihood ----
# Parameterization: b_ps, b_sp free; a_p, a_s fixed structural.
# All n-1 bivariate Gaussian transition densities share the same F, Q
# (dt constant), so these are computed once per call.

ou_coupled_log_lik <- function(b_ps, b_sp, sigma_p, sigma_s, c_p, c_s,
                               x_mat, u_vec, a_p, a_s, dt,
                               n_q_steps = 40L) {
  if (a_p * a_s - b_ps * b_sp <= 0) return(-Inf)   # stability gate

  mats <- tryCatch(
    ou_coupled_matrices(a_p, a_s, b_ps, b_sp, sigma_p, sigma_s, dt, n_q_steps),
    error = function(e) NULL
  )
  if (is.null(mats) || !is.finite(mats$log_det_Q)) return(-Inf)

  n      <- nrow(x_mat) - 1L
  c_vec  <- c(c_p, c_s)
  ll_base <- -n * log(2 * pi) - (n / 2) * mats$log_det_Q

  quad <- 0
  for (i in seq_len(n)) {
    mu_i <- as.vector(mats$F %*% x_mat[i, ])
    u_i  <- u_vec[i]
    if (u_i != 0 && !is.null(mats$A_inv))
      mu_i <- mu_i + as.vector((mats$F - diag(2L)) %*% mats$A_inv %*% (c_vec * u_i))
    r    <- x_mat[i + 1L, ] - mu_i
    quad <- quad + as.numeric(t(r) %*% mats$Q_inv %*% r)
  }
  ll_base - 0.5 * quad
}

# ---- Numerical MLE for the coupled OU block ----
# Optimises (b_ps, b_sp, log σ_p, log σ_s, c_p, c_s).
# Initialises from the closed-form uncoupled estimates for σ and c.

ou_coupled_mle <- function(x_mat, u_vec, a_p, a_s, dt,
                           b_ps_init = 0, b_sp_init = 0) {
  init_p <- ou_mle(x_mat[, 1L], u_vec, a_p, dt)
  init_s <- ou_mle(x_mat[, 2L], u_vec, a_s, dt)

  th0 <- c(b_ps_init, b_sp_init,
           log(max(init_p$sigma, 1e-5)), log(max(init_s$sigma, 1e-5)),
           init_p$c_gain, init_s$c_gain)

  neg_ll <- function(th) {
    ll <- tryCatch(
      ou_coupled_log_lik(th[1L], th[2L], exp(th[3L]), exp(th[4L]),
                         th[5L], th[6L],
                         x_mat, u_vec, a_p, a_s, dt),
      error = function(e) -Inf)
    if (is.finite(ll)) -ll else 1e10
  }

  res <- optim(th0, neg_ll, method = "L-BFGS-B",
               control = list(maxit = 300L, factr = 1e8))
  th  <- res$par
  sp_hat <- exp(th[3L]);  ss_hat <- exp(th[4L])

  # Numerical Hessian → standard errors
  eps <- 1e-4;  np <- length(th)
  H   <- matrix(0, np, np)
  for (i in seq_len(np)) for (j in i:np) {
    ei <- ej <- numeric(np); ei[i] <- eps; ej[j] <- eps
    H[i, j] <- (neg_ll(th+ei+ej) - neg_ll(th+ei-ej) -
                  neg_ll(th-ei+ej) + neg_ll(th-ei-ej)) / (4 * eps^2)
    H[j, i] <- H[i, j]
  }
  V       <- tryCatch(solve(H), error = function(e) matrix(NA_real_, np, np))
  safe_se <- function(k) sqrt(max(V[k, k], 0, na.rm = TRUE))

  list(
    b_ps = th[1L], b_sp = th[2L],
    sigma_p = sp_hat, sigma_s = ss_hat,
    c_p = th[5L], c_s = th[6L],
    b_ps_se    = safe_se(1L), b_sp_se    = safe_se(2L),
    sigma_p_se = sp_hat * safe_se(3L),
    sigma_s_se = ss_hat * safe_se(4L),
    c_p_se     = safe_se(5L), c_s_se     = safe_se(6L),
    log_lik    = -res$value,
    convergence = res$convergence,
    n_obs = nrow(x_mat) - 1L
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
  w       <- exp(delta_vec)          # w_k = 1/g_k

  # Closed-form MLE for mu_0
  denom_w <- sum(w)
  if (!is.finite(denom_w) || denom_w < .Machine$double.eps) {
    warning("ig_obs_mle: denominator sum(exp(delta)) near zero; returning NAs.")
    return(list(mu0=NA_real_, kappa=NA_real_, mu0_se=NA_real_, kappa_se=NA_real_,
                log_lik=NA_real_, n_beats=length(tau_vec)))
  }
  mu0_hat <- sum(tau_vec * w^2) / denom_w

  # Closed-form MLE for kappa (given mu_0_hat)
  mu_k    <- mu0_hat * g
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
  i_log_mu0   <- kappa_hat * denom_w / mu0_hat   # reuse denom_w from above
  i_log_kappa <- n / 2
  # SEs are unreliable if kappa is at its cap; flag with NA
  mu0_se   <- if (kappa_hat >= 1e6 - 1) NA_real_ else mu0_hat   / sqrt(i_log_mu0)
  kappa_se <- if (kappa_hat >= 1e6 - 1) NA_real_ else kappa_hat / sqrt(i_log_kappa)

  ll <- sum(log_ig_pdf(tau_vec, mu_k, kappa_hat))

  list(mu0 = mu0_hat, kappa = max(kappa_hat, 1e-4),
       mu0_se = mu0_se, kappa_se = kappa_se,
       log_lik = ll, n_beats = n)
}

# ---- Combined conditional MLE for all free parameters ----

full_conditional_mle <- function(sim_res) {
  sp    <- sim_res$params$structural
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  n_t   <- length(sim_res$time)
  u_vec <- extract_u_vec(sim_res)

  fp          <- sim_res$params$free
  use_coupled <- (abs(fp$b_ps) + abs(fp$b_sp)) > 1e-12

  if (use_coupled) {
    x_mat  <- cbind(sim_res$p, sim_res$s)
    mle_ou <- ou_coupled_mle(x_mat, u_vec, sp$a_p, sp$a_s, dt,
                             b_ps_init = fp$b_ps, b_sp_init = fp$b_sp)
    # Wrap into the same interface expected by the rest of the function
    mle_p  <- list(sigma    = mle_ou$sigma_p, c_gain   = mle_ou$c_p,
                   sigma_se = mle_ou$sigma_p_se, c_se  = mle_ou$c_p_se,
                   log_lik  = NA_real_)
    mle_s  <- list(sigma    = mle_ou$sigma_s, c_gain   = mle_ou$c_s,
                   sigma_se = mle_ou$sigma_s_se, c_se  = mle_ou$c_s_se,
                   log_lik  = NA_real_)
  } else {
    mle_ou <- NULL
    mle_p  <- ou_mle(sim_res$p, u_vec, sp$a_p, dt)
    mle_s  <- ou_mle(sim_res$s, u_vec, sp$a_s, dt)
  }

  spk <- sim_res$spikes

  if (length(spk) < 2L) {
    na6 <- setNames(rep(NA_real_, 6),
                    c("sigma_p","sigma_s","mu0","kappa","c_p","c_s"))
    return(c(as.list(na6), list(n_beats = 0L)))
  }
  tau_vec  <- diff(spk)

  # Compute mean delta over the actual interval instead of the endpoint
  delta_v <- compute_effective_delta(spk, sim_res$time, sim_res$delta)

  mle_obs  <- ig_obs_mle(tau_vec, delta_v)

  ll_ou <- if (!is.null(mle_ou)) {
    mle_ou$log_lik
  } else {
    mle_p$log_lik + mle_s$log_lik
  }

  list(
    b_ps     = if (!is.null(mle_ou)) mle_ou$b_ps    else fp$b_ps,
    b_sp     = if (!is.null(mle_ou)) mle_ou$b_sp    else fp$b_sp,
    b_ps_se  = if (!is.null(mle_ou)) mle_ou$b_ps_se else NA_real_,
    b_sp_se  = if (!is.null(mle_ou)) mle_ou$b_sp_se else NA_real_,
    sigma_p  = mle_p$sigma,   sigma_p_se  = mle_p$sigma_se,
    c_p      = mle_p$c_gain,  c_p_se      = mle_p$c_se,
    sigma_s  = mle_s$sigma,   sigma_s_se  = mle_s$sigma_se,
    c_s      = mle_s$c_gain,  c_s_se      = mle_s$c_se,
    mu0      = mle_obs$mu0,   mu0_se      = mle_obs$mu0_se,
    kappa    = mle_obs$kappa, kappa_se    = mle_obs$kappa_se,
    n_beats  = mle_obs$n_beats,
    ll_obs   = mle_obs$log_lik,
    ll_total = ll_ou + mle_obs$log_lik
  )
}

# ---- Observed information matrix (negative Hessian) for IG block ----
#
# Parameterization: theta = (log mu_0, log kappa) for numerical stability.
# Returns the 2x2 FIM, its inverse (parameter covariance), SE, and
# condition number.

ig_obs_fim <- function(mu0, kappa_val, tau_vec, delta_vec, eps = 1e-5) {
  th0 <- c(log(mu0), log(kappa_val))

  obj <- function(lth) {
    mk <- exp(lth[1L]) * exp(-delta_vec)
    if (any(mk <= 0)) return(-Inf)
    sum(log_ig_pdf(tau_vec, mk, exp(lth[2L])))
  }

  # Symmetric numerical Hessian
  H <- matrix(0, 2, 2)
  for (i in 1:2) {
    for (j in i:2) {
      ei <- ej <- rep(0, 2)
      ei[i] <- eps; ej[j] <- eps
      H[i, j] <- (obj(th0+ei+ej) - obj(th0+ei-ej) -
                  obj(th0-ei+ej) + obj(th0-ei-ej)) / (4 * eps^2)
      H[j, i] <- H[i, j]
    }
  }
  fim <- -H
  fim_inv <- tryCatch(solve(fim), error = function(e) matrix(NA, 2, 2))
  list(
    fim       = fim,
    se        = sqrt(abs(diag(fim_inv))),  # SE in log-space
    cond      = tryCatch(kappa(fim), error = function(e) NA_real_),
    eig       = tryCatch(eigen(fim, symmetric = TRUE, only.values = TRUE)$values,
                         error = function(e) rep(NA_real_, 2))
  )
}

# ---- Observed information matrix for OU block ----
#
# Parameterization: theta = (log sigma, c_gain).
# The cross-term is exactly zero at the MLE (by the score equation),
# so the matrix is diagonal.

ou_fim <- function(sigma_val, x_vec, u_vec, a, dt) {
  n    <- length(x_vec) - 1L
  z    <- u_vec[seq_len(n)] / a * (-expm1(-a * dt))
  C_a  <- (-expm1(-2 * a * dt)) / (2 * a)
  ssz  <- sum(z^2)

  # I(log sigma) = 2*n  [information = 2 per observation in log-sigma scale]
  # I(c_gain)    = sum(z^2) / (sigma^2 * C_a)
  i_log_sig <- 2 * n
  i_c       <- if (ssz > 1e-14) ssz / (sigma_val^2 * C_a) else 0

  fim <- diag(c(i_log_sig, i_c))
  list(
    fim       = fim,
    se        = c(1 / sqrt(i_log_sig), if (i_c > 0) 1 / sqrt(i_c) else Inf),
    cond      = if (i_c > 0) max(i_log_sig, i_c) / min(i_log_sig, i_c) else Inf,
    eig       = c(i_log_sig, i_c)
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

.full_fim_coupled <- function(sim_res, mle, sp, u_vec, dt) {
  x_mat   <- cbind(sim_res$p, sim_res$s)
  spk     <- sim_res$spikes
  tau_vec <- diff(spk)
  delta_v <- compute_effective_delta(spk, sim_res$time, sim_res$delta)

  th0 <- c(mle$b_ps, mle$b_sp,
           log(mle$sigma_p), log(mle$sigma_s),
           mle$c_p, mle$c_s,
           log(mle$mu0), log(mle$kappa))

  obj <- function(th) {
    if (sp$a_p * sp$a_s - th[1L] * th[2L] <= 0) return(-Inf)
    ll_ou <- tryCatch(
      ou_coupled_log_lik(th[1L], th[2L], exp(th[3L]), exp(th[4L]),
                         th[5L], th[6L],
                         x_mat, u_vec, sp$a_p, sp$a_s, dt),
      error = function(e) -Inf)
    mu_k   <- exp(th[7L]) * exp(-delta_v)
    ll_obs <- sum(log_ig_pdf(tau_vec, mu_k, exp(th[8L])))
    ll_ou + ll_obs
  }

  eps <- 1e-4;  np <- length(th0)
  H   <- matrix(0, np, np)
  for (i in seq_len(np)) for (j in i:np) {
    ei <- ej <- numeric(np); ei[i] <- eps; ej[j] <- eps
    H[i, j] <- (obj(th0+ei+ej) - obj(th0+ei-ej) -
                  obj(th0-ei+ej) + obj(th0-ei-ej)) / (4 * eps^2)
    H[j, i] <- H[i, j]
  }
  FIM <- -H
  pnames <- c("b_ps", "b_sp",
              "log(sigma_p)", "log(sigma_s)", "c_p", "c_s",
              "log(mu_0)", "log(kappa)")
  dimnames(FIM) <- list(pnames, pnames)
  eig_all <- tryCatch(
    eigen(FIM, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, np))
  list(FIM         = FIM,
       eigenvalues = eig_all,
       cond_number = max(eig_all) / min(eig_all[eig_all > 0]),
       param_names = pnames,
       blocks      = list(p   = list(eig = NA, cond = NA),
                          s   = list(eig = NA, cond = NA),
                          obs = list(eig = tail(eig_all, 2L), cond = NA)))
}

full_conditional_fim <- function(sim_res) {
  sp    <- sim_res$params$structural
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  n_t   <- length(sim_res$time)
  u_vec <- extract_u_vec(sim_res)

  mle <- full_conditional_mle(sim_res)

  if ((abs(sim_res$params$free$b_ps) + abs(sim_res$params$free$b_sp)) > 1e-12)
    return(.full_fim_coupled(sim_res, mle, sp, u_vec, dt))

  fim_p   <- ou_fim(mle$sigma_p, sim_res$p, u_vec, sp$a_p, dt)
  fim_s   <- ou_fim(mle$sigma_s, sim_res$s, u_vec, sp$a_s, dt)

  spk      <- sim_res$spikes
  tau_vec  <- diff(spk)
  delta_v <- compute_effective_delta(spk, sim_res$time, sim_res$delta)
  fim_obs  <- ig_obs_fim(mle$mu0, mle$kappa, tau_vec, delta_v)

  # Assemble 6x6 block-diagonal matrix
  FIM <- matrix(0, 6, 6)
  FIM[1:2, 1:2] <- fim_p$fim
  FIM[3:4, 3:4] <- fim_s$fim
  FIM[5:6, 5:6] <- fim_obs$fim

  param_names <- c("log(sigma_p)", "c_p",
                   "log(sigma_s)", "c_s",
                   "log(mu_0)",    "log(kappa)")
  dimnames(FIM) <- list(param_names, param_names)

  eig_all <- tryCatch(
    eigen(FIM, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, 6)
  )

  list(
    FIM         = FIM,
    eigenvalues = eig_all,
    cond_number = max(eig_all) / min(eig_all[eig_all > 0]),
    se          = c(fim_p$se, fim_s$se, fim_obs$se),
    param_names = param_names,
    blocks      = list(p = fim_p, s = fim_s, obs = fim_obs)
  )
}

# ---- Profile likelihood for each free parameter ----
#
# For each parameter, the profile is computed by maximizing over the
# remaining parameters in its block. Thanks to block separability,
# this reduces to a 1-D optimization in each case.

ou_log_lik_at <- function(x_vec, u_vec, a, sigma, c_gain, dt) {
  n       <- length(x_vec) - 1L
  idx     <- seq_len(n)
  e_adt   <- exp(-a * dt)
  means   <- x_vec[idx] * e_adt + (c_gain * u_vec[idx] / a) * (-expm1(-a * dt))
  var_v   <- (sigma^2 / (2 * a)) * (-expm1(-2 * a * dt))
  sum(dnorm(x_vec[idx + 1L], mean = means, sd = sqrt(var_v), log = TRUE))
}

profile_lik_one <- function(param, grid, sim_res, mle = NULL) {
  sp    <- sim_res$params$structural
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  n_t   <- length(sim_res$time)
  u_vec <- extract_u_vec(sim_res)

  spk     <- sim_res$spikes
  tau_vec <- diff(spk)
  delta_v <- compute_effective_delta(spk, sim_res$time, sim_res$delta)

  # --- Pre-compute nuisance MLEs once (they do not depend on the grid value) ---
  # OU blocks: c_hat does not depend on sigma; sigma_hat DOES depend on c.
  # IG block:  mu0_hat does not depend on kappa; kappa_hat depends on mu0.

  # Parasympathetic block
  n_p  <- length(sim_res$p) - 1L
  e_p  <- exp(-sp$a_p * dt)
  z_p  <- u_vec[seq_len(n_p)] / sp$a_p * (-expm1(-sp$a_p * dt))
  y_p  <- sim_res$p[-1L] - sim_res$p[seq_len(n_p)] * e_p
  C_ap <- (-expm1(-2 * sp$a_p * dt)) / (2 * sp$a_p)
  c_p_hat <- if (sum(z_p^2) > 1e-14) sum(y_p * z_p) / sum(z_p^2) else 0
  # sigma_hat at joint MLE (used for sigma_p profile; c_hat is nuisance there)
  r_p_joint  <- y_p - c_p_hat * z_p
  sig_p_hat  <- sqrt(max(mean(r_p_joint^2) / C_ap, 1e-14))

  # Sympathetic block
  n_s  <- length(sim_res$s) - 1L
  e_s  <- exp(-sp$a_s * dt)
  z_s  <- u_vec[seq_len(n_s)] / sp$a_s * (-expm1(-sp$a_s * dt))
  y_s  <- sim_res$s[-1L] - sim_res$s[seq_len(n_s)] * e_s
  C_as <- (-expm1(-2 * sp$a_s * dt)) / (2 * sp$a_s)
  c_s_hat <- if (sum(z_s^2) > 1e-14) sum(y_s * z_s) / sum(z_s^2) else 0
  r_s_joint  <- y_s - c_s_hat * z_s
  sig_s_hat  <- sqrt(max(mean(r_s_joint^2) / C_as, 1e-14))

  # IG block: mu0_hat is independent of kappa
  g       <- exp(-delta_v)
  w       <- exp(delta_v)
  mu0_hat <- sum(tau_vec * w^2) / sum(w)
  mu_k_at_mu0hat <- mu0_hat * g

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
           kappa = {
             # Profile mu0 over fixed kappa = v (mu0_hat is kappa-independent)
             if (v <= 0) return(-Inf)
             sum(log_ig_pdf(tau_vec, mu_k_at_mu0hat, v))
           },
           sigma_p = {
             # Profile c_p over fixed sigma_p = v (c_hat is sigma-independent)
             if (v <= 0) return(-Inf)
             ou_log_lik_at(sim_res$p, u_vec, sp$a_p, v, c_p_hat, dt)
           },
           sigma_s = {
             if (v <= 0) return(-Inf)
             ou_log_lik_at(sim_res$s, u_vec, sp$a_s, v, c_s_hat, dt)
           },
           c_p = {
             # Profile sigma_p over fixed c_p = v (sigma_hat DOES depend on v)
             r_v   <- y_p - v * z_p
             sig_v <- sqrt(max(mean(r_v^2) / C_ap, 1e-14))
             ou_log_lik_at(sim_res$p, u_vec, sp$a_p, sig_v, v, dt)
           },
           c_s = {
             r_v   <- y_s - v * z_s
             sig_v <- sqrt(max(mean(r_v^2) / C_as, 1e-14))
             ou_log_lik_at(sim_res$s, u_vec, sp$a_s, sig_v, v, dt)
           },
           b_ps = {
             # Fix b_ps = v; profile over (b_sp, log σ_p, log σ_s, c_p, c_s).
             if (sp$a_p * sp$a_s - v * sim_res$params$free$b_sp <= 0) return(-Inf)
             x_mat <- cbind(sim_res$p, sim_res$s)
             mle_v <- if (!is.null(mle)) mle else full_conditional_mle(sim_res)
             th0   <- c(mle_v$b_sp,
                        log(max(mle_v$sigma_p, 1e-5)),
                        log(max(mle_v$sigma_s, 1e-5)),
                        mle_v$c_p, mle_v$c_s)
             nll <- function(th) {
               if (sp$a_p * sp$a_s - v * th[1L] <= 0) return(1e10)
               ll <- tryCatch(
                 ou_coupled_log_lik(v, th[1L], exp(th[2L]), exp(th[3L]),
                                    th[4L], th[5L],
                                    x_mat, u_vec, sp$a_p, sp$a_s, dt),
                 error = function(e) -Inf)
               if (is.finite(ll)) -ll else 1e10
             }
             -optim(th0, nll, method = "L-BFGS-B",
                    control = list(maxit = 150L, factr = 1e9))$value
           },
           b_sp = {
             if (sp$a_p * sp$a_s - sim_res$params$free$b_ps * v <= 0) return(-Inf)
             x_mat <- cbind(sim_res$p, sim_res$s)
             mle_v <- if (!is.null(mle)) mle else full_conditional_mle(sim_res)
             th0   <- c(mle_v$b_ps,
                        log(max(mle_v$sigma_p, 1e-5)),
                        log(max(mle_v$sigma_s, 1e-5)),
                        mle_v$c_p, mle_v$c_s)
             nll <- function(th) {
               if (sp$a_p * sp$a_s - th[1L] * v <= 0) return(1e10)
               ll <- tryCatch(
                 ou_coupled_log_lik(th[1L], v, exp(th[2L]), exp(th[3L]),
                                    th[4L], th[5L],
                                    x_mat, u_vec, sp$a_p, sp$a_s, dt),
                 error = function(e) -Inf)
               if (is.finite(ll)) -ll else 1e10
             }
             -optim(th0, nll, method = "L-BFGS-B",
                    control = list(maxit = 150L, factr = 1e9))$value
           },
           -Inf   # unknown parameter name
    )
  }, numeric(1L))
}

# ---- All profile likelihoods in one call ----

all_profile_likelihoods <- function(sim_res, n_grid = 60, width = 3.0) {
  mle <- full_conditional_mle(sim_res)

  fp_sim      <- sim_res$params$free
  use_coupled <- (abs(fp_sim$b_ps) + abs(fp_sim$b_sp)) > 1e-12

  specs <- list(
    mu0     = list(center = mle$mu0,     log_scale = TRUE,
                   label = expression(mu[0] ~ "(s)")),
    kappa   = list(center = mle$kappa,   log_scale = TRUE,
                   label = expression(kappa)),
    sigma_p = list(center = mle$sigma_p, log_scale = TRUE,
                   label = expression(sigma[p])),
    sigma_s = list(center = mle$sigma_s, log_scale = TRUE,
                   label = expression(sigma[s])),
    c_p     = list(center = mle$c_p,     log_scale = FALSE, width_override = 5.0,
                   label = expression(c[p])),
    c_s     = list(center = mle$c_s,     log_scale = FALSE, width_override = 5.0,
                   label = expression(c[s]))
  )
  if (use_coupled) {
    specs$b_ps <- list(center = mle$b_ps, log_scale = FALSE, width_override = 3.0,
                       label = expression(b[ps]))
    specs$b_sp <- list(center = mle$b_sp, log_scale = FALSE, width_override = 3.0,
                       label = expression(b[sp]))
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
    ll <- profile_lik_one(p, grid, sim_res, mle = mle)
    list(param = p, grid = grid, ll = ll,
         mle = s$center, label = s$label,
         ll_mle = max(ll, na.rm = TRUE))
  }) |> setNames(names(specs))
}
