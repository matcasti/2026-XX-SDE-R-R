# ============================================================
# R/mle.R
#
# MLE for the SDE-IG model.
# Requires: model_functions.R sourced first.
#
# Log-likelihood separability still holds:
#
#   L(theta | states) = L_p(sigma_p, c_p)    [OU parasympathetic — exact closed form]
#                     + L_s(sigma_s, c_s)    [OU sympathetic    — exact closed form]
#                     + L_obs(mu_0, kappa)   [IG observations   — exact PP likelihood]
#
# For the OU blocks the transition density is exactly Gaussian,
# giving closed-form MLEs with no approximation error.
#
# For the IG observation block, treating inter-beat intervals as i.i.d.
# IG(mu_0 * exp(-delta_k), kappa) is biased whenever mu(t) varies within
# an interval — which it always does for the fast vagal OU branch
# (time constant tau_p = 0.5 s vs typical RR interval ~0.85 s, giving
# within-interval autocorrelation e^{-2*0.85} ~ 0.18).
#
# The exact conditional log-likelihood for the doubly-stochastic point process is
#
#   L_obs(mu_0, kappa) = sum_k  log lambda(t_k)
#                      - int_0^T lambda(t) dt
#
# where lambda(t) = ig_hazard(tau(t), mu_0 * exp(-delta(t)), kappa) and
# tau(t) = elapsed time since the most recent spike before t.
# This is maximised numerically (L-BFGS-B, 2 parameters) using the
# full stored trajectory. The closed-form density-product MLEs are
# retained as ig_obs_mle_approx() and used only as warm-start values.
# ============================================================

# ---- Block 1 / 2: Exact closed-form OU MLE ----

ou_mle <- function(x_vec, u_vec, a, dt) {
  n     <- length(x_vec) - 1L
  if (n < 2L) return(list(sigma = NA_real_, c_gain = NA_real_,
                          sigma_se = NA_real_, c_se = NA_real_,
                          log_lik = NA_real_))
  e_adt <- exp(-a * dt)
  z     <- u_vec[seq_len(n)] / a * (-expm1(-a * dt))
  y     <- x_vec[-1L] - x_vec[seq_len(n)] * e_adt

  ssz   <- sum(z^2)
  c_hat <- if (ssz > 1e-14) sum(y * z) / ssz else 0
  resid <- y - c_hat * z
  C_a       <- (-expm1(-2 * a * dt)) / (2 * a)
  sigma2    <- max(mean(resid^2) / C_a, 1e-14)
  sigma_hat <- sqrt(sigma2)
  sigma_se  <- sigma_hat / sqrt(2 * n)
  c_se      <- if (ssz > 1e-14) sigma_hat * sqrt(C_a / ssz) else Inf
  ll        <- sum(dnorm(resid, 0, sigma_hat * sqrt(C_a), log = TRUE))

  list(sigma = sigma_hat, c_gain = c_hat,
       sigma_se = sigma_se, c_se = c_se,
       log_lik = ll, n_obs = n)
}

# ---- Vectorised OU log-likelihood (for OU profiles) ----

ou_log_lik_at <- function(x_vec, u_vec, a, sigma, c_gain, dt) {
  n     <- length(x_vec) - 1L
  e_adt <- exp(-a * dt)
  C_a   <- (-expm1(-2 * a * dt)) / (2 * a)
  z     <- u_vec[seq_len(n)] / a * (-expm1(-a * dt))
  means <- x_vec[seq_len(n)] * e_adt + c_gain * z
  sum(dnorm(x_vec[-1L], means, sigma * sqrt(C_a), log = TRUE))
}

# ---- IG density-based approximate MLE (warm start only) ----
#
# Uses START-of-interval delta (less biased than end-of-interval).
# Closed-form score equations; biased when mu(t) varies within intervals.

ig_obs_mle_approx <- function(tau_vec, delta_vec) {
  stopifnot(length(tau_vec) == length(delta_vec), length(tau_vec) > 1L)
  g         <- exp(-delta_vec)
  w         <- exp(delta_vec)
  mu0_hat   <- sum(tau_vec * w^2) / sum(w)
  mu_k      <- mu0_hat * g
  psi       <- sum((tau_vec - mu_k)^2 / (mu_k^2 * tau_vec))
  n         <- length(tau_vec)
  kappa_hat <- n / max(psi, 1e-10)
  list(mu0 = mu0_hat, kappa = max(kappa_hat, 1e-4))
}

# ---- Exact conditional PP log-likelihood for IG block ----
#
# L_obs(log_mu0, log_kappa) = sum_k log lambda(t_k)  [event term]
#                           - dt * sum_i lambda(t_i)  [integral term]
#
# Arguments:
#   log_theta  : c(log mu_0, log kappa)
#   elapsed    : tau(t_i) at every grid point
#   g_all      : exp(-delta(t_i))  [= mu(t_i) / mu_0]
#   dt         : step size
#   beat_mask  : logical, TRUE at grid points that are beat times
#   valid_mask : logical, TRUE where elapsed > 1e-6 (exclude resets)

ig_obs_pp_log_lik <- function(log_theta,
                              elapsed, g_all, dt,
                              beat_mask, valid_mask) {
  mu0   <- exp(log_theta[1L])
  kappa <- exp(log_theta[2L])

  # Integral term
  mu_v   <- mu0 * g_all[valid_mask]
  el_v   <- elapsed[valid_mask]
  lam_v  <- exp(log_ig_hazard_v(el_v, mu_v, kappa))
  lam_v[!is.finite(lam_v)] <- 0
  ll_int <- -dt * sum(lam_v)

  # Event term
  bv     <- beat_mask & valid_mask
  ll_ev  <- sum(log_ig_hazard_v(elapsed[bv], mu0 * g_all[bv], kappa))

  if (!is.finite(ll_ev) || !is.finite(ll_int)) return(-Inf)
  ll_ev + ll_int
}

# ---- Exact conditional MLE for IG block ----
#
# Optimises ig_obs_pp_log_lik over (log mu_0, log kappa).
# Returns point estimates, delta-method SEs, and the 2x2 observed FIM.

ig_obs_mle <- function(sim_res) {
  tg    <- sim_res$time
  dt    <- tg[2L] - tg[1L]
  spk   <- sim_res$spikes
  delta <- sim_res$delta

  if (length(spk) < 2L)
    return(list(mu0 = NA_real_, kappa = NA_real_,
                mu0_se = NA_real_, kappa_se = NA_real_,
                log_lik = NA_real_, n_beats = 0L,
                converged = FALSE, fim_obs = matrix(NA_real_, 2, 2)))

  # Precompute (done once per call)
  elapsed    <- compute_elapsed_times(tg, spk)
  g_all      <- exp(-delta)
  beat_idx   <- vapply(spk, function(t) which.min(abs(tg - t)), integer(1L))
  beat_mask  <- logical(length(tg)); beat_mask[beat_idx] <- TRUE
  valid_mask <- elapsed > 1e-6

  # Warm start: approx MLE with start-of-interval delta
  tau_vec  <- diff(spk)
  si_idx   <- vapply(spk[-length(spk)],
                     function(t) which.min(abs(tg - t)), integer(1L))
  init     <- ig_obs_mle_approx(tau_vec, delta[si_idx])
  start    <- c(log(max(init$mu0,   1e-4)),
                log(max(init$kappa, 0.1 )))

  neg_ll <- function(lth) {
    ll <- tryCatch(
      ig_obs_pp_log_lik(lth, elapsed, g_all, dt, beat_mask, valid_mask),
      error = function(e) -Inf)
    if (is.finite(ll)) -ll else 1e10
  }

  opt <- tryCatch(
    optim(start, neg_ll, method = "L-BFGS-B",
          lower = c(-5, -2), upper = c(5, 7),
          control = list(maxit = 500L, factr = 1e7)),
    error = function(e)
      list(par = start, value = Inf, convergence = 99L))

  mu0_hat   <- exp(opt$par[1L])
  kappa_hat <- exp(opt$par[2L])

  # Observed FIM: negative Hessian of log-lik via symmetric finite differences
  eps <- 1e-4
  H   <- matrix(0, 2, 2)
  for (ii in 1:2) for (jj in ii:2) {
    ei <- ej <- c(0, 0); ei[ii] <- eps; ej[jj] <- eps
    lth <- opt$par
    pp_ll <- function(lt)
      ig_obs_pp_log_lik(lt, elapsed, g_all, dt, beat_mask, valid_mask)
    H[ii, jj] <- (pp_ll(lth+ei+ej) - pp_ll(lth+ei-ej) -
                    pp_ll(lth-ei+ej) + pp_ll(lth-ei-ej)) / (4 * eps^2)
    H[jj, ii] <- H[ii, jj]
  }
  fim_obs <- -H
  fim_inv <- tryCatch(solve(fim_obs), error = function(e) matrix(NA_real_, 2, 2))
  se_log  <- sqrt(abs(diag(fim_inv)))

  list(
    mu0       = mu0_hat,
    kappa     = kappa_hat,
    mu0_se    = mu0_hat   * se_log[1L],
    kappa_se  = kappa_hat * se_log[2L],
    log_lik   = -opt$value,
    n_beats   = length(spk),
    converged = (opt$convergence == 0L),
    fim_obs   = fim_obs
  )
}

# ---- Combined conditional MLE (all free parameters) ----

full_conditional_mle <- function(sim_res) {
  sp    <- sim_res$params$structural
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  n_t   <- length(sim_res$time)
  u_vec <- if (!is.null(sim_res$input_fn))
    vapply(sim_res$time, sim_res$input_fn, numeric(1L))
  else rep(0, n_t)

  mle_p   <- ou_mle(sim_res$p, u_vec, sp$a_p, dt)
  mle_s   <- ou_mle(sim_res$s, u_vec, sp$a_s, dt)
  mle_obs <- ig_obs_mle(sim_res)

  if (is.na(mle_obs$mu0)) {
    na6 <- setNames(rep(NA_real_, 6),
                    c("sigma_p","sigma_s","mu0","kappa","c_p","c_s"))
    return(c(as.list(na6),
             list(n_beats = 0L, fim_obs = matrix(NA_real_, 2, 2))))
  }

  list(
    sigma_p  = mle_p$sigma,    sigma_p_se = mle_p$sigma_se,
    c_p      = mle_p$c_gain,   c_p_se     = mle_p$c_se,
    sigma_s  = mle_s$sigma,    sigma_s_se = mle_s$sigma_se,
    c_s      = mle_s$c_gain,   c_s_se     = mle_s$c_se,
    mu0      = mle_obs$mu0,    mu0_se     = mle_obs$mu0_se,
    kappa    = mle_obs$kappa,  kappa_se   = mle_obs$kappa_se,
    n_beats  = mle_obs$n_beats,
    ll_p     = mle_p$log_lik,
    ll_s     = mle_s$log_lik,
    ll_obs   = mle_obs$log_lik,
    ll_total = mle_p$log_lik + mle_s$log_lik + mle_obs$log_lik,
    fim_obs  = mle_obs$fim_obs   # 2x2 exact PP FIM for (log mu_0, log kappa)
  )
}

# ---- Observed information matrix for OU block ----

ou_fim <- function(sigma_val, c_gain_val, x_vec, u_vec, a, dt) {
  n   <- length(x_vec) - 1L
  z   <- u_vec[seq_len(n)] / a * (-expm1(-a * dt))
  C_a <- (-expm1(-2 * a * dt)) / (2 * a)
  ssz <- sum(z^2)
  i_log_sig <- 2 * n
  i_c       <- if (ssz > 1e-14) ssz / (sigma_val^2 * C_a) else 0
  fim <- diag(c(i_log_sig, i_c))
  list(fim  = fim,
       se   = c(1 / sqrt(i_log_sig), if (i_c > 0) 1 / sqrt(i_c) else Inf),
       cond = if (i_c > 0) i_log_sig / i_c else Inf,
       eig  = c(i_log_sig, i_c))
}

# ---- Full 6x6 block-diagonal FIM ----
#
# IG block: uses exact PP Hessian returned by ig_obs_mle (via full_conditional_mle).
# OU blocks: analytic diagonal FIM.

full_conditional_fim <- function(sim_res) {
  sp    <- sim_res$params$structural
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  n_t   <- length(sim_res$time)
  u_vec <- if (!is.null(sim_res$input_fn))
    vapply(sim_res$time, sim_res$input_fn, numeric(1L))
  else rep(0, n_t)

  mle     <- full_conditional_mle(sim_res)
  fim_p   <- ou_fim(mle$sigma_p, mle$c_p, sim_res$p, u_vec, sp$a_p, dt)
  fim_s   <- ou_fim(mle$sigma_s, mle$c_s, sim_res$s, u_vec, sp$a_s, dt)
  fim_obs_blk <- mle$fim_obs

  fim_obs_inv <- tryCatch(solve(fim_obs_blk),
                          error = function(e) matrix(NA_real_, 2, 2))
  obs_se  <- sqrt(abs(diag(fim_obs_inv)))
  obs_eig <- tryCatch(
    eigen(fim_obs_blk, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, 2))
  obs_cond <- if (all(is.finite(obs_eig)) && min(obs_eig) > 0)
    max(obs_eig) / min(obs_eig) else NA_real_

  FIM <- matrix(0, 6, 6)
  FIM[1:2, 1:2] <- fim_p$fim
  FIM[3:4, 3:4] <- fim_s$fim
  FIM[5:6, 5:6] <- fim_obs_blk

  param_names <- c("log(sigma_p)", "c_p",
                   "log(sigma_s)", "c_s",
                   "log(mu_0)",    "log(kappa)")
  dimnames(FIM) <- list(param_names, param_names)

  eig_all <- tryCatch(
    eigen(FIM, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, 6))

  list(FIM         = FIM,
       eigenvalues = eig_all,
       cond_number = max(eig_all) / min(eig_all[eig_all > 0]),
       se          = c(fim_p$se, fim_s$se, obs_se),
       param_names = param_names,
       blocks      = list(p   = fim_p,
                          s   = fim_s,
                          obs = list(fim  = fim_obs_blk,
                                     se   = obs_se,
                                     cond = obs_cond,
                                     eig  = obs_eig)))
}

# ---- Profile likelihood: OU parameters (exact analytic) ----

ou_profile_one <- function(param, grid, sim_res) {
  sp    <- sim_res$params$structural
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  n_t   <- length(sim_res$time)
  u_vec <- if (!is.null(sim_res$input_fn))
    vapply(sim_res$time, sim_res$input_fn, numeric(1L))
  else rep(0, n_t)

  vapply(grid, function(v) {
    if (!is.finite(v)) return(-Inf)
    switch(param,
           sigma_p = { mle_c <- ou_mle(sim_res$p, u_vec, sp$a_p, dt)
           ou_log_lik_at(sim_res$p, u_vec, sp$a_p, v, mle_c$c_gain, dt) },
           sigma_s = { mle_c <- ou_mle(sim_res$s, u_vec, sp$a_s, dt)
           ou_log_lik_at(sim_res$s, u_vec, sp$a_s, v, mle_c$c_gain, dt) },
           c_p     = { mle_s <- ou_mle(sim_res$p, u_vec, sp$a_p, dt)
           ou_log_lik_at(sim_res$p, u_vec, sp$a_p, mle_s$sigma, v, dt) },
           c_s     = { mle_s <- ou_mle(sim_res$s, u_vec, sp$a_s, dt)
           ou_log_lik_at(sim_res$s, u_vec, sp$a_s, mle_s$sigma, v, dt) },
           -Inf)
  }, numeric(1L))
}

# ---- Profile likelihood: IG parameters (quadratic from exact PP FIM) ----
#
# For mu_0 and kappa, the profile log-likelihood is constructed from the
# observed FIM at the exact PP MLE via the quadratic approximation:
#
#   L_profile(theta_i; theta_hat) = -0.5 * (log theta_i - log theta_hat_i)^2
#                                    * I_partial_ii
#
# where I_partial_ii = FIM[i,i] - FIM[i,j]^2 / FIM[j,j] is the profile
# (partial) information for parameter i after profiling out j.
# The FIM is nearly diagonal for this model, so I_partial ~ FIM[i,i].
# The quadratic profile is appropriate here because the recovery study
# confirms asymptotic normality of the exact PP MLE (coverage ~ 95%).

ig_quadratic_profile <- function(param, grid_orig, mle_val, fim_obs) {
  i        <- if (param == "mu0") 1L else 2L
  j        <- 3L - i
  I_ii     <- fim_obs[i, i]
  I_ij     <- fim_obs[i, j]
  I_jj     <- fim_obs[j, j]
  I_partial <- I_ii - if (abs(I_jj) > 1e-10) I_ij^2 / I_jj else 0
  log_hat  <- log(mle_val)
  -0.5 * (log(grid_orig) - log_hat)^2 * I_partial
}

# ---- All profile likelihoods in one call ----

all_profile_likelihoods <- function(sim_res, n_grid = 60, width = 3.0) {
  mle     <- full_conditional_mle(sim_res)
  fim_obs <- mle$fim_obs

  make_grid <- function(center, log_scale) {
    if (log_scale)
      exp(seq(log(center / width), log(center * width), length.out = n_grid))
    else {
      spread <- max(abs(center) * width, 0.5)
      seq(center - spread, center + spread, length.out = n_grid)
    }
  }

  # OU parameters: exact analytic profiles
  ou_specs <- list(
    sigma_p = list(center = mle$sigma_p, log_scale = TRUE,
                   label = expression(sigma[p])),
    sigma_s = list(center = mle$sigma_s, log_scale = TRUE,
                   label = expression(sigma[s])),
    c_p     = list(center = mle$c_p,     log_scale = FALSE,
                   label = expression(c[p])),
    c_s     = list(center = mle$c_s,     log_scale = FALSE,
                   label = expression(c[s]))
  )
  ou_pf <- lapply(names(ou_specs), function(p) {
    s    <- ou_specs[[p]]
    grid <- make_grid(s$center, s$log_scale)
    ll   <- ou_profile_one(p, grid, sim_res)
    ll_n <- ll - max(ll, na.rm = TRUE)
    list(param = p, grid = grid, ll = ll_n,
         mle = s$center, label = s$label, ll_mle = 0)
  }) |> setNames(names(ou_specs))

  # IG parameters: quadratic profile from exact PP FIM
  ig_specs <- list(
    mu0   = list(center = mle$mu0,   label = expression(mu[0] ~ "(s)")),
    kappa = list(center = mle$kappa, label = expression(kappa))
  )
  ig_pf <- lapply(names(ig_specs), function(p) {
    s    <- ig_specs[[p]]
    grid <- make_grid(s$center, TRUE)
    ll_n <- ig_quadratic_profile(p, grid, s$center, fim_obs)
    list(param = p, grid = grid, ll = ll_n,
         mle = s$center, label = s$label, ll_mle = 0)
  }) |> setNames(names(ig_specs))

  # Canonical order expected by plot_profiles()
  c(ig_pf["mu0"], ig_pf["kappa"],
    ou_pf["sigma_p"], ou_pf["sigma_s"],
    ou_pf["c_p"],     ou_pf["c_s"])
}
