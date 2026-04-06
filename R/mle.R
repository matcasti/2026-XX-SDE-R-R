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
    if (length(idx) > 0) {
      -log(mean(exp(-delta_vec[idx])))
    } else {
      delta_vec[which.min(abs(time_vec - t_end))]
    }
  }, numeric(1L))
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

  g       <- exp(-delta_vec)         # g_k = mu_k / mu_0  (=1 when delta=0)
  w       <- exp(delta_vec)          # w_k = 1/g_k

  # Closed-form MLE for mu_0
  mu0_hat <- sum(tau_vec * w^2) / sum(w)

  # Closed-form MLE for kappa (given mu_0_hat)
  mu_k    <- mu0_hat * g
  psi     <- sum((tau_vec - mu_k)^2 / (mu_k^2 * tau_vec))
  n       <- length(tau_vec)
  kappa_hat <- n / max(psi, 1e-10)

  # Standard errors via analytic FIM (derived in supplementary):
  #   I(log mu_0) = kappa * sum(exp(delta_k)) / mu_0
  #   I(log kappa) = N / 2
  #   Cross term = 0 (block diagonal in log-space)
  i_log_mu0   <- kappa_hat * sum(w) / mu0_hat
  i_log_kappa <- n / 2

  mu0_se   <- mu0_hat   / sqrt(i_log_mu0)
  kappa_se <- kappa_hat / sqrt(i_log_kappa)

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
  u_vec <- if (!is.null(sim_res$input_fn))
    vapply(sim_res$time, sim_res$input_fn, numeric(1L))
  else rep(0, n_t)

  mle_p <- ou_mle(sim_res$p, u_vec, sp$a_p, dt)
  mle_s <- ou_mle(sim_res$s, u_vec, sp$a_s, dt)

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

  list(
    sigma_p  = mle_p$sigma,   sigma_p_se  = mle_p$sigma_se,
    c_p      = mle_p$c_gain,  c_p_se      = mle_p$c_se,
    sigma_s  = mle_s$sigma,   sigma_s_se  = mle_s$sigma_se,
    c_s      = mle_s$c_gain,  c_s_se      = mle_s$c_se,
    mu0      = mle_obs$mu0,   mu0_se      = mle_obs$mu0_se,
    kappa    = mle_obs$kappa, kappa_se    = mle_obs$kappa_se,
    n_beats  = mle_obs$n_beats,
    ll_p     = mle_p$log_lik,
    ll_s     = mle_s$log_lik,
    ll_obs   = mle_obs$log_lik,
    ll_total = mle_p$log_lik + mle_s$log_lik + mle_obs$log_lik
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

full_conditional_fim <- function(sim_res) {
  sp    <- sim_res$params$structural
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  n_t   <- length(sim_res$time)
  u_vec <- if (!is.null(sim_res$input_fn))
    vapply(sim_res$time, sim_res$input_fn, numeric(1L))
  else rep(0, n_t)

  mle <- full_conditional_mle(sim_res)

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
  n <- length(x_vec) - 1L
  sum(vapply(seq_len(n), function(i)
    ou_log_transition_density(x_vec[i+1L], x_vec[i], a, sigma, c_gain, u_vec[i], dt),
    numeric(1L)))
}

profile_lik_one <- function(param, grid, sim_res) {
  sp    <- sim_res$params$structural
  dt    <- sim_res$time[2L] - sim_res$time[1L]
  n_t   <- length(sim_res$time)
  u_vec <- if (!is.null(sim_res$input_fn))
    vapply(sim_res$time, sim_res$input_fn, numeric(1L))
  else rep(0, n_t)

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
             # Profile kappa over fixed mu0 = v (closed form)
             if (v <= 0) return(-Inf)
             mu_k <- v * g
             psi  <- sum((tau_vec - mu_k)^2 / (mu_k^2 * tau_vec))
             k    <- length(tau_vec) / max(psi, 1e-10)
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
           -Inf   # unknown parameter name
    )
  }, numeric(1L))
}

# ---- All profile likelihoods in one call ----

all_profile_likelihoods <- function(sim_res, n_grid = 60, width = 3.0) {
  mle <- full_conditional_mle(sim_res)

  specs <- list(
    mu0     = list(center = mle$mu0,     log_scale = TRUE,
                   label = expression(mu[0] ~ "(s)")),
    kappa   = list(center = mle$kappa,   log_scale = TRUE,
                   label = expression(kappa)),
    sigma_p = list(center = mle$sigma_p, log_scale = TRUE,
                   label = expression(sigma[p])),
    sigma_s = list(center = mle$sigma_s, log_scale = TRUE,
                   label = expression(sigma[s])),
    c_p     = list(center = mle$c_p,     log_scale = FALSE,
                   label = expression(c[p])),
    c_s     = list(center = mle$c_s,     log_scale = FALSE,
                   label = expression(c[s]))
  )

  lapply(names(specs), function(p) {
    s <- specs[[p]]
    grid <- if (s$log_scale) {
      exp(seq(log(s$center / width), log(s$center * width), length.out = n_grid))
    } else {
      spread <- max(abs(s$center) * width, 0.5)
      seq(s$center - spread, s$center + spread, length.out = n_grid)
    }
    ll <- profile_lik_one(p, grid, sim_res)
    list(param = p, grid = grid, ll = ll,
         mle = s$center, label = s$label,
         ll_mle = max(ll, na.rm = TRUE))
  }) |> setNames(names(specs))
}
