# ============================================================
# R/api.R
#
# High-level user-facing API for the SDE-IG framework.
#
# Quick start
# -----------
#   source every R/*.R file once, then:
#
#   fit <- sde_ig(rr_vec)          # rr_vec: numeric vector of IBIs in seconds
#   print(fit)                      # parameter table + fit statistics
#   autoplot(fit)                   # tachogram / Δ̂(t) / innovations (3 panels)
#   autoplot(fit, "states")         # Δ̂(t) + p̂(t), ŝ(t) branches
#   autoplot(fit, "spectrum")       # theoretical PSD with LF/HF bands
#   diagnose(fit)                   # time-rescaling 4-panel diagnostics
#   get_states(fit)                 # data.frame of filtered state trajectory
#   augment(fit)                    # per-beat data.frame (.fitted, .resid, …)
#
# Soft dependencies for plotting: ggplot2 (>= 3.4), patchwork.
# Hard dependencies: all other R/*.R files sourced first.
# ============================================================

# ---- 0. Internal helpers --------------------------------------------------------

.check_gg <- function() {
  if (!requireNamespace("ggplot2",   quietly = TRUE))
    stop("ggplot2 required for plots: install.packages('ggplot2')")
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("patchwork required for multi-panel plots: install.packages('patchwork')")
  invisible(TRUE)
}

.theme_sde_ig <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.35),
      axis.line        = ggplot2::element_line(colour = "grey70", linewidth = 0.35),
      strip.text       = ggplot2::element_text(face = "bold", size = 10),
      plot.title       = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.35, "cm"),
      legend.text      = ggplot2::element_text(size = 9)
    )
}

.COL <- c(
  delta  = "#2C5F2E",   # net drive Δ
  p_hat  = "#0072B2",   # parasympathetic
  s_hat  = "#D55E00",   # sympathetic
  mu     = "#1a1a1a",   # fitted mean interval
  rr     = "grey60",    # observed R-R
  innov  = "#5B5EA6",   # innovations
  warn   = "#D55E00"    # KS / CI bands
)

# ---- 1. Constructor -------------------------------------------------------------

#' Fit the SDE-IG autonomic model to R-R interval data
#'
#' @param rr        Numeric vector of **inter-beat intervals in seconds**
#'                  (positive, finite). Alternatively pass cumulative spike
#'                  times and set `type = "spikes"`.
#' @param type      `"rr"` (default) or `"spikes"`.
#' @param coupled   Logical. Also estimate cross-inhibition terms
#'                  \eqn{a_{ps}} and \eqn{a_{sp}}? Default `FALSE`.
#' @param input_fn  Optional function `u(t)` describing a known exogenous
#'                  protocol (e.g. tilt table, exercise ramp). If `NULL`
#'                  (default) the input is treated as identically zero.
#' @param params_init Optional `make_model_params()` output with starting
#'                  values. `NULL` triggers automatic spectral initialisation.
#' @param method    `"twostage"` (default) or `"mle"`.  The two-stage
#'                  wavelet-Whittle + concentrated-MLE pipeline is more
#'                  robust to spectral ridges and is recommended.
#' @param verbose   Print L-BFGS-B progress? Default `FALSE`.
#'
#' @return An `"sde_ig_fit"` object.  Key slots:
#' \describe{
#'   \item{`free_hat`}{Named list of parameter estimates.}
#'   \item{`ll`, `aic`, `bic`}{Model fit statistics.}
#'   \item{`filter`}{UKF run at the MLE (passed to `autoplot` / `diagnose`).}
#'   \item{`convergence`}{L-BFGS-B code (0 = clean convergence).}
#' }
#'
#' @examples
#' \dontrun{
#' ## Simulate 5-min recording and fit
#' p   <- make_model_params(a_p=2, a_s=0.2, sigma_p=0.20, sigma_s=0.07,
#'                          mu_0=1.0, rho=0.20)
#' sim <- sim_sde_ig(300, 0.005, p)
#' fit <- sde_ig(diff(sim$spikes))
#' print(fit)
#' autoplot(fit)
#' diagnose(fit)
#'
#' ## Real data: load RR intervals from a plain-text file (one per line, seconds)
#' rr  <- scan("subject01_rr.txt")
#' fit <- sde_ig(rr)
#' coef(fit)
#' augment(fit)       # per-beat data frame for downstream analysis
#' }
sde_ig <- function(rr,
                   type        = c("rr", "spikes"),
                   coupled     = FALSE,
                   input_fn    = NULL,
                   params_init = NULL,
                   method      = c("twostage", "mle"),
                   verbose     = FALSE) {

  mc     <- match.call()
  type   <- match.arg(type)
  method <- match.arg(method)

  # ---- normalise input ----
  if (type == "spikes") {
    sp <- as.numeric(rr)
    if (!all(diff(sp) > 0))
      stop("sde_ig: spike times must be strictly increasing.")
    rr_vec <- diff(sp)
    spikes <- sp
  } else {
    rr_vec <- as.numeric(rr)
    bad <- !is.finite(rr_vec) | rr_vec <= 0
    if (any(bad))
      stop(sprintf("sde_ig: %d non-positive or non-finite values in `rr`.",
                   sum(bad)))
    if (length(rr_vec) < 30L)
      warning("sde_ig: fewer than 30 beats — estimates may be unreliable.")
    spikes <- c(0, cumsum(rr_vec))
  }

  if (is.null(input_fn)) input_fn <- function(t) 0

  # ---- starting values ----
  if (is.null(params_init)) {
    params_init <- tryCatch(
      spectral_init(rr_vec, verbose = verbose),
      error = function(e) {
        if (verbose)
          message("sde_ig: spectral_init failed (", conditionMessage(e),
                  "); using generic defaults.")
        make_model_params()
      }
    )
  }
  if (coupled && (abs(params_init$free$a_ps) + abs(params_init$free$a_sp)) < 1e-12) {
    cpl <- tryCatch(
      band_filtered_coupling_init(rr_vec,
                                  a_p = params_init$free$a_p,
                                  a_s = params_init$free$a_s),
      error = function(e) list(a_ps = 0.05, a_sp = 0.05)
    )
    fp0 <- params_init$free
    params_init <- make_model_params(
      a_p     = fp0$a_p,   a_s     = fp0$a_s,
      a_ps    = max(cpl$a_ps, 1e-4),
      a_sp    = max(cpl$a_sp, 1e-4),
      sigma_p = fp0$sigma_p, sigma_s = fp0$sigma_s,
      mu_0    = fp0$mu_0,  rho     = fp0$rho
    )
  }

  # ---- fit ----
  estimator <- if (method == "twostage") pp_mle_twostage else pp_mle

  result <- tryCatch(
    estimator(spikes, params_init, input_fn, verbose = verbose),
    error = function(e)
      stop("sde_ig: optimisation failed — ", conditionMessage(e),
           "\nTry method = 'mle', or supply explicit params_init.")
  )

  if (is.null(result) || is.na(result$ll))
    warning("sde_ig: log-likelihood is NA — check fit$convergence and ",
            "consider different starting values.")

  # ---- information criteria (6 or 8 natural parameters) ----
  np  <- if (coupled) 8L else 6L
  ll  <- result$ll
  aic <- if (is.finite(ll)) -2 * ll + 2 * np           else NA_real_
  bic <- if (is.finite(ll)) -2 * ll + log(length(rr_vec)) * np else NA_real_

  structure(
    list(
      call        = mc,
      rr          = rr_vec,
      spikes      = spikes,
      n_beats     = length(rr_vec),
      params_hat  = result$params_hat,
      free_hat    = result$free_hat,
      ll          = ll,
      aic         = aic,
      bic         = bic,
      convergence = result$convergence,
      filter      = result$filter,
      method      = method,
      coupled     = coupled,
      input_fn    = input_fn,
      se_hat      = NULL          # populated lazily by sde_ig_add_se()
    ),
    class = "sde_ig_fit"
  )
}

# ---- 2. Standard S3 generics ----------------------------------------------------

#' @export
print.sde_ig_fit <- function(x, digits = 4L, ...) {
  if (!is.numeric(digits)) {
    stop(sprintf(
      paste0("print.sde_ig_fit: `digits` must be numeric, not \"%s\".\n",
             "  For visualization use: plot(fit, type = \"%s\") ",
             "or autoplot(fit, \"%s\")"),
      digits, digits, digits), call. = FALSE)
  }
  conv_msg <- if (x$convergence == 0L) "ok (code 0)"
  else sprintf("WARNING (code %d) — re-run or check starting values",
               x$convergence)
  cat("\n\u2500\u2500 SDE-IG Autonomic Model",
      "\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500",
      "\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  Method : %-12s  Coupled : %-5s  Convergence : %s\n",
              x$method, if (x$coupled) "yes" else "no", conv_msg))
  cat(sprintf("  Beats  : %-6d        Duration : %.1f s\n",
              x$n_beats, tail(x$spikes, 1L)))
  cat(sprintf("  Log-lik: %10.2f   AIC : %10.2f   BIC : %.2f\n\n",
              x$ll, x$aic, x$bic))
  cat("  Parameter estimates\n")
  print(.coef_table(x), digits = digits, row.names = FALSE, right = FALSE)
  if (is.null(x$se_hat))
    cat("\n  Standard errors: call sde_ig_add_se(fit) [slow]\n")
  invisible(x)
}

.coef_table <- function(x) {
  fp    <- x$free_hat
  kap   <- kappa_from_rho(fp$mu_0, fp$rho)
  pnms  <- c("a_p [Hz]", "a_s [Hz]", "sigma_p", "sigma_s",
             "mu_0 [s]", "rho", "kappa")
  vals  <- c(fp$a_p, fp$a_s, fp$sigma_p, fp$sigma_s, fp$mu_0, fp$rho, kap)
  interp <- c(
    sprintf("Vagal decay      tau_p = %.2f s", 1 / fp$a_p),
    sprintf("Sympath. decay   tau_s = %.1f s", 1 / fp$a_s),
    "Vagal noise amplitude",
    "Sympathetic noise amplitude",
    sprintf("Baseline interval  HR = %.1f bpm", 60 / fp$mu_0),
    "Baseline coefficient of variation",
    "IG shape  (derived: mu_0/rho^2)"
  )
  if (!is.null(x$se_hat)) {
    se  <- x$se_hat
    se_v <- c(se$a_p_se, se$a_s_se, se$sigma_p_se, se$sigma_s_se,
              NA_real_, se$rho_se, NA_real_)
    df <- data.frame(Parameter = pnms, Estimate = round(vals, 5),
                     SE = round(se_v, 5), Interpretation = interp,
                     stringsAsFactors = FALSE)
  } else {
    df <- data.frame(Parameter = pnms, Estimate = round(vals, 5),
                     Interpretation = interp, stringsAsFactors = FALSE)
  }
  if (x$coupled) {
    extra_vals  <- c(fp$a_ps, fp$a_sp)
    extra_interp <- c("s \u2192 p inhibition", "p \u2192 s inhibition")
    extra_se    <- if (!is.null(x$se_hat))
      c(x$se_hat$a_ps_se, x$se_hat$a_sp_se) else NULL
    extra <- if (!is.null(extra_se))
      data.frame(Parameter = c("a_ps", "a_sp"),
                 Estimate = round(extra_vals, 5),
                 SE = round(extra_se, 5),
                 Interpretation = extra_interp, stringsAsFactors = FALSE)
    else
      data.frame(Parameter = c("a_ps", "a_sp"),
                 Estimate = round(extra_vals, 5),
                 Interpretation = extra_interp, stringsAsFactors = FALSE)
    df <- rbind(df, extra)
  }
  df
}

#' @export
summary.sde_ig_fit <- function(object, ...) {
  print(object, ...)
  fp   <- object$free_hat
  var_p <- fp$sigma_p^2 / (2 * fp$a_p)
  var_s <- fp$sigma_s^2 / (2 * fp$a_s)
  cat("\n  Derived spectral quantities\n")
  cat(sprintf("    Vagal stationary var   (A_p): %.5f s^2\n", var_p))
  cat(sprintf("    Sympathetic stat. var  (A_s): %.5f s^2\n", var_s))
  cat(sprintf("    LF/HF power ratio          : %.3f\n",
              (fp$sigma_s^2 / fp$a_s) / (fp$sigma_p^2 / fp$a_p)))
  cat(sprintf("    HF pole at freq            : %.4f Hz\n", fp$a_p / (2 * pi)))
  cat(sprintf("    LF pole at freq            : %.4f Hz\n", fp$a_s / (2 * pi)))
  invisible(object)
}

#' @export
coef.sde_ig_fit <- function(object, ...) {
  fp <- object$free_hat
  v  <- c(a_p = fp$a_p, a_s = fp$a_s,
          sigma_p = fp$sigma_p, sigma_s = fp$sigma_s,
          mu_0 = fp$mu_0, rho = fp$rho)
  if (object$coupled) v <- c(v, a_ps = fp$a_ps, a_sp = fp$a_sp)
  v
}

#' @export
logLik.sde_ig_fit <- function(object, ...) {
  np <- if (object$coupled) 8L else 6L
  structure(object$ll, df = np, nobs = object$n_beats, class = "logLik")
}

#' @export
AIC.sde_ig_fit <- function(object, ..., k = 2) {
  -2 * object$ll + k * (if (object$coupled) 8L else 6L)
}

#' @export
BIC.sde_ig_fit <- function(object, ...) {
  -2 * object$ll + log(object$n_beats) * (if (object$coupled) 8L else 6L)
}

# ---- 3. Broom-compatible generics -----------------------------------------------

#' Tidy parameter table (broom-compatible)
#' @export
tidy.sde_ig_fit <- function(x, ...) {
  df    <- .coef_table(x)
  out   <- df[, intersect(c("Parameter", "Estimate", "SE"), names(df))]
  names(out)[1L] <- "term"
  names(out)[2L] <- "estimate"
  out
}

#' One-row model summary (broom-compatible)
#' @export
glance.sde_ig_fit <- function(x, ...) {
  data.frame(
    n_beats     = x$n_beats,
    duration_s  = round(tail(x$spikes, 1L), 2),
    logLik      = round(x$ll,  2),
    AIC         = round(x$aic, 2),
    BIC         = round(x$bic, 2),
    convergence = x$convergence,
    method      = x$method,
    coupled     = x$coupled,
    stringsAsFactors = FALSE
  )
}

#' Per-beat augmented data frame (broom-compatible)
#'
#' Returns one row per inter-beat interval with the original IBI, the
#' filter-predicted mean interval, the raw innovation, and the posterior
#' state estimates at each beat time.
#' @export
augment.sde_ig_fit <- function(x, ...) {
  flt <- x$filter
  data.frame(
    beat_time  = flt$beat_times,
    rr         = x$rr,
    .fitted    = flt$mu_filt,
    .resid     = flt$innov,
    .std_resid = flt$innov / sqrt(flt$S_innov),
    delta_hat  = flt$delta_filt,
    p_hat      = flt$m_filt[, "p"],
    s_hat      = flt$m_filt[, "s"],
    stringsAsFactors = FALSE
  )
}

# ---- 4. Standard errors and confidence intervals --------------------------------

#' Add standard errors to an sde_ig_fit object (slow)
#'
#' Computes a numerical Hessian of the marginal log-likelihood at the MLE and
#' stores the resulting per-parameter SEs in `fit$se_hat`.  Expect 15–90 s
#' for a 300-beat recording (one UKF run per Hessian entry).
#'
#' SEs are approximate: they derive from a first-order delta method mapping
#' optimizer coordinates `(log a_s, log(a_p-a_s), log A_p, log A_s, log kappa)`
#' back to natural parameters.  Use `all_profile_likelihoods()` for more
#' accurate intervals on the conditional likelihood.
#'
#' @param fit  An `sde_ig_fit` object.
#' @param eps  Step size for numerical differentiation (default `1e-4`).
#' @return The same object with `se_hat` populated.
sde_ig_add_se <- function(fit, eps = 1e-4) {
  stopifnot(inherits(fit, "sde_ig_fit"))
  fp       <- fit$free_hat
  inp      <- fit$input_fn
  mu_bar   <- mean(fit$rr)
  coupled  <- fit$coupled

  # Mirror pp_mle pack() at the MLE to obtain the optimizer-coordinate vector.
  A_p   <- max(fp$sigma_p^2 / (2 * fp$a_p), 1e-12)
  A_s   <- max(fp$sigma_s^2 / (2 * fp$a_s), 1e-12)
  ratio <- max(fp$a_p / max(fp$a_s, 1e-8), SDE_IG_MIN_AP_AS_RATIO)
  mu_0_anchor <- max(mu_bar * exp(-(A_p + A_s) / 2), 0.10)
  kap   <- kappa_from_rho(fp$mu_0, fp$rho)
  th    <- c(log(fp$a_s), log(ratio), log(A_p), log(A_s), log(max(kap, 0.5)))
  if (coupled) {
    c_v <- max(fp$a_ps * fp$a_sp, 1e-8)
    r_v <- max(fp$a_ps / max(fp$a_sp, 1e-8), 1e-3)
    th  <- c(th, log(c_v), log(r_v))
  }

  # Local unpack mirrors pp_mle's closure.
  # NOTE: must be updated if pp_mle's parameterisation changes.
  .unp <- function(v) {
    a_s_v  <- exp(v[1L]); a_p_v <- a_s_v + exp(v[2L])
    A_p_v  <- max(exp(v[3L]), 1e-10); A_s_v <- max(exp(v[4L]), 1e-10)
    kap_v  <- max(exp(v[5L]), 0.5)
    mu0_v  <- max(mu_bar * exp(-(A_p_v + A_s_v) / 2), 0.10)
    a_ps_v <- if (coupled) sqrt(max(exp(v[6L]) * exp(v[7L]), 0)) else 0
    a_sp_v <- if (coupled) sqrt(max(exp(v[6L]) / exp(v[7L]), 0)) else 0
    list(a_p = a_p_v, a_s = a_s_v, a_ps = a_ps_v, a_sp = a_sp_v,
         sigma_p = sqrt(max(2 * a_p_v * A_p_v, 1e-10)),
         sigma_s = sqrt(max(2 * a_s_v * A_s_v, 1e-10)),
         mu_0 = mu0_v, rho = sqrt(mu0_v / kap_v),
         c_p = fp$c_p, c_s = fp$c_s)
  }

  neg_ll_v <- function(v) {
    params_v <- list(structural = fit$params_hat$structural, free = .unp(v))
    flt <- tryCatch(pp_ukf(fit$spikes, params_v, inp),
                    error = function(e) list(ll_vec = NULL))
    if (is.null(flt$ll_vec)) return(1e10)
    lf <- flt$ll_vec[is.finite(flt$ll_vec)]
    if (length(lf) < 2L) return(1e10)
    -sum(lf)
  }

  np <- length(th)
  message(sprintf(
    "sde_ig_add_se: computing %d\u00d7%d Hessian (%d UKF calls)  …",
    np, np, 2L * np * (np + 1L) / 2L))
  H <- numerical_hessian(neg_ll_v, th, eps = eps)
  V <- tryCatch(solve(H), error = function(e) {
    warning("sde_ig_add_se: Hessian singular — SEs are NA.")
    matrix(NA_real_, np, np)
  })
  se_v <- sqrt(pmax(diag(V), 0, na.rm = TRUE))   # SE on log / natural scale

  # Delta-method back-transform to natural-scale SEs (first-order).
  fit$se_hat <- list(
    a_s_se     = fp$a_s  * se_v[1L],   # d(a_s)/d(log a_s) = a_s
    a_p_se     = fp$a_p  * se_v[2L],   # d(a_p)/d(log ratio) = a_p  [first-order]
    sigma_p_se = fp$sigma_p / 2     * se_v[3L],   # sigma_p ~ A_p^0.5
    sigma_s_se = fp$sigma_s / 2     * se_v[4L],
    rho_se     = fp$rho     / 2     * se_v[5L],   # rho ~ kappa^-0.5
    a_ps_se    = if (coupled) fp$a_ps * se_v[6L] else NA_real_,
    a_sp_se    = if (coupled) fp$a_sp * se_v[7L] else NA_real_
  )
  message("Done.")
  fit
}

#' @export
confint.sde_ig_fit <- function(object, parm, level = 0.95,
                               method = c("hessian", "none"), ...) {
  method <- match.arg(method)
  if (method == "none") {
    message("confint: use method = 'hessian' to compute approximate Wald CIs.")
    return(invisible(NULL))
  }
  if (is.null(object$se_hat)) object <- sde_ig_add_se(object)

  z    <- qnorm((1 + level) / 2)
  fp   <- object$free_hat
  se   <- object$se_hat
  ests <- c(fp$a_p,     fp$a_s,     fp$sigma_p, fp$sigma_s, fp$rho)
  ses  <- c(se$a_p_se, se$a_s_se, se$sigma_p_se, se$sigma_s_se, se$rho_se)
  pnms <- c("a_p", "a_s", "sigma_p", "sigma_s", "rho")
  if (object$coupled) {
    ests <- c(ests, fp$a_ps, fp$a_sp)
    ses  <- c(ses,  se$a_ps_se, se$a_sp_se)
    pnms <- c(pnms, "a_ps", "a_sp")
  }
  pct <- function(p) sprintf("%.1f%%", 100 * p)
  out <- data.frame(
    term     = pnms,
    estimate = ests,
    lo       = pmax(ests - z * ses, 0),
    hi       = ests + z * ses,
    row.names = NULL, stringsAsFactors = FALSE
  )
  names(out)[3:4] <- c(pct((1 - level) / 2), pct((1 + level) / 2))
  out
}

# ---- 5. Data accessors ----------------------------------------------------------

#' Extract filtered state trajectory as a tidy data frame
#'
#' @param fit  An `sde_ig_fit` object.
#' @param dt   Grid resolution in seconds (default 0.005).
#' @return A data.frame with columns `time`, `delta_hat`, `sd_delta`,
#'   `p_hat`, `s_hat`, `mu_hat`.
get_states <- function(fit, dt = 0.005) {
  stopifnot(inherits(fit, "sde_ig_fit"))
  tg  <- seq(0, tail(fit$spikes, 1L), by = dt)
  fg  <- filter_to_grid(fit$filter, tg)
  data.frame(time     = tg,
             delta_hat = fg$delta,
             sd_delta  = fg$sd_delta,
             p_hat     = fg$p,
             s_hat     = fg$s,
             mu_hat    = fg$mu)
}

#' Extract per-beat tachogram with fitted mean interval
#'
#' @param fit  An `sde_ig_fit` object.
#' @return A data.frame with columns `beat_time`, `rr`, `mu_hat`, `delta_hat`.
get_tachogram <- function(fit) {
  stopifnot(inherits(fit, "sde_ig_fit"))
  flt <- fit$filter
  data.frame(beat_time = flt$beat_times,
             rr        = fit$rr,
             mu_hat    = flt$mu_filt,
             delta_hat = flt$delta_filt)
}

# ---- 6. Internal plot builders --------------------------------------------------

.build_tachogram <- function(fit) {
  td <- get_tachogram(fit)
  ggplot2::ggplot(td, ggplot2::aes(x = beat_time)) +
    ggplot2::geom_point(ggplot2::aes(y = rr),
                        colour = .COL["rr"], size = 0.55, alpha = 0.55) +
    ggplot2::geom_line(ggplot2::aes(y = mu_hat),
                       colour = .COL["mu"], linewidth = 0.85) +
    ggplot2::labs(x = "Time (s)", y = "R\u2013R interval (s)",
                  title = "Tachogram",
                  subtitle = expression(paste("Dots: observed  |  line: fitted ",
                                              hat(mu)(t)))) +
    .theme_sde_ig()
}

.build_delta <- function(sdf) {
  ggplot2::ggplot(sdf, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = delta_hat - 2 * sd_delta,
                                      ymax = delta_hat + 2 * sd_delta),
                         fill = .COL["delta"], alpha = 0.15, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = delta_hat),
                       colour = .COL["delta"], linewidth = 0.85, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.35) +
    ggplot2::labs(x = NULL,
                  y = expression(hat(Delta)(t) == hat(s)(t) - hat(p)(t)),
                  title = "Net autonomic drive",
                  subtitle = "\u00b12 posterior SD band") +
    .theme_sde_ig()
}

.build_branches <- function(sdf) {
  n    <- nrow(sdf)
  dl   <- rbind(
    data.frame(time = sdf$time, value = sdf$p_hat,
               branch = "Parasympathetic \u2014 p\u0302(t)", stringsAsFactors = FALSE),
    data.frame(time = sdf$time, value = sdf$s_hat,
               branch = "Sympathetic \u2014 s\u0302(t)",     stringsAsFactors = FALSE)
  )
  lvls <- c("Parasympathetic \u2014 p\u0302(t)", "Sympathetic \u2014 s\u0302(t)")
  dl$branch <- factor(dl$branch, levels = lvls)
  cols <- stats::setNames(c(.COL["p_hat"], .COL["s_hat"]), lvls)

  ggplot2::ggplot(dl, ggplot2::aes(x = time, y = value, colour = branch)) +
    ggplot2::geom_line(linewidth = 0.75, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.35) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "Time (s)", y = "Latent state",
                  title = "Branch state estimates",
                  subtitle = "UKF posterior means") +
    .theme_sde_ig()
}

.build_innovations <- function(fit) {
  flt <- fit$filter
  idf <- data.frame(time  = flt$beat_times,
                    innov = flt$innov,
                    sd_s  = sqrt(flt$S_innov))
  ggplot2::ggplot(idf, ggplot2::aes(x = time)) +
    ggplot2::geom_segment(ggplot2::aes(xend = time,
                                       y = -2 * sd_s, yend = 2 * sd_s),
                          colour = .COL["innov"], alpha = 0.25, linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(y = innov),
                        colour = .COL["innov"], size = 0.65, alpha = 0.55) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.35) +
    ggplot2::labs(x = "Time (s)",
                  y = expression(tau[k] - hat(mu)[k]),
                  title = "Filter innovations",
                  subtitle = "Points: residuals  |  bars: \u00b12\u221aS\u2096") +
    .theme_sde_ig()
}

.build_spectrum <- function(fit) {
  fp    <- fit$free_hat
  freqs <- seq(0.005, 0.50, length.out = 600)
  lor   <- function(a, sig) sig^2 / (a^2 + (2 * pi * freqs)^2)
  psd_p <- lor(fp$a_p, fp$sigma_p)
  psd_s <- lor(fp$a_s, fp$sigma_s)
  lvls  <- c("Total", "Vagal (HF)", "Sympathetic (LF)")
  df    <- data.frame(
    freq      = rep(freqs, 3L),
    psd       = c(psd_p + psd_s, psd_p, psd_s),
    component = factor(rep(lvls, each = 600L), levels = lvls)
  )
  cols <- stats::setNames(c("#1a1a1a", .COL["p_hat"], .COL["s_hat"]), lvls)
  ltys <- stats::setNames(c("solid", "dashed", "dotdash"), lvls)
  ymax <- max(psd_p + psd_s)

  ggplot2::ggplot(df, ggplot2::aes(x = freq, y = psd,
                                   colour = component, linetype = component)) +
    ggplot2::annotate("rect", xmin = 0.04, xmax = 0.15,
                      ymin = -Inf, ymax = Inf,
                      alpha = 0.05, fill = .COL["s_hat"]) +
    ggplot2::annotate("rect", xmin = 0.15, xmax = 0.40,
                      ymin = -Inf, ymax = Inf,
                      alpha = 0.05, fill = .COL["p_hat"]) +
    ggplot2::annotate("text", x = 0.095, y = ymax * 0.88,
                      label = "LF", colour = .COL["s_hat"],
                      size = 3.2, fontface = "bold") +
    ggplot2::annotate("text", x = 0.275,  y = ymax * 0.88,
                      label = "HF", colour = .COL["p_hat"],
                      size = 3.2, fontface = "bold") +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_linetype_manual(values = ltys, name = NULL) +
    ggplot2::scale_x_continuous(
      breaks = c(0.04, 0.15, 0.40),
      labels = c("0.04", "0.15", "0.40")) +
    ggplot2::labs(
      x = "Frequency (Hz)",
      y = expression(paste("PSD of ",Delta(t),"  [s"^2,"/Hz]")),
      title = "Theoretical power spectral density",
      subtitle = sprintf(
        "Lorentzian poles at %.3f Hz (HF) and %.4f Hz (LF)",
        fp$a_p / (2 * pi), fp$a_s / (2 * pi))) +
    .theme_sde_ig()
}

# ---- 7. plot / autoplot ---------------------------------------------------------

#' Plot an sde_ig_fit object
#'
#' Wrapper around `autoplot`; returns a patchwork object.
#'
#' @param x    An `sde_ig_fit` object.
#' @param type One of `"overview"`, `"states"`, `"tachogram"`, `"spectrum"`.
#' @param dt   State-grid resolution in seconds (default 0.005).
#' @param ...  Forwarded to `autoplot`.
#' @export
plot.sde_ig_fit <- function(x, type = "overview", dt = 0.005, ...) {
  .check_gg()
  autoplot(x, type = type, dt = dt, ...)
}

#' ggplot2 autoplot for sde_ig_fit
#'
#' @param object An `sde_ig_fit` object.
#' @param type   Plot type:
#' \describe{
#'   \item{`"overview"`}{(default) Tachogram + Δ̂(t) ± SD + innovations (3 rows).}
#'   \item{`"states"`}{Δ̂(t) + branch estimates p̂(t), ŝ(t) (2 rows).}
#'   \item{`"tachogram"`}{Observed RR with fitted μ(t) (single panel).}
#'   \item{`"spectrum"`}{Theoretical PSD with LF/HF bands (single panel).}
#' }
#' @param dt     State-grid resolution in seconds (default 0.005).
#' @param ...    Ignored.
#' @return A `ggplot` / `patchwork` object; pipe-modify with `& ggplot2::theme(...)`.
#' @export
autoplot.sde_ig_fit <- function(object,
                                type = c("overview", "states",
                                         "tachogram", "spectrum"),
                                dt = 0.005, ...) {
  .check_gg()
  type <- match.arg(type)
  switch(type,
         overview = {
           sdf <- get_states(object, dt)
           p1  <- .build_tachogram(object)
           p2  <- .build_delta(sdf)
           p3  <- .build_innovations(object)
           patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(1.2, 1, 0.9))
         },
         states = {
           sdf <- get_states(object, dt)
           p1  <- .build_delta(sdf) + ggplot2::labs(x = NULL)
           p2  <- .build_branches(sdf)
           patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 1.1))
         },
         tachogram = .build_tachogram(object),
         spectrum  = .build_spectrum(object)
  )
}

# ---- 8. Time-rescaling diagnostics ----------------------------------------------

#' S3 generic for model diagnostics
#' @export
diagnose <- function(x, ...) UseMethod("diagnose")

#' Time-rescaling diagnostics for sde_ig_fit
#'
#' Implements the Brown et al. (2002) four-panel suite:
#' (1) Uniform Q-Q plot of \eqn{u_k = 1-\exp(-\Lambda_k)};
#' (2) Normal-score Q-Q;
#' (3) ACF of \eqn{u_k};
#' (4) Empirical CDF of \eqn{\Lambda_k} vs Exp(1).
#'
#' @param x    An `sde_ig_fit` object.
#' @param ...  Ignored.
#' @return A patchwork object (4 panels, 2 × 2).
#' @export
diagnose.sde_ig_fit <- function(x, ...) {
  .check_gg()

  # Build sim_res proxy for compute_time_rescaling
  dt  <- 0.005
  tg  <- seq(0, tail(x$spikes, 1L), by = dt)
  fg  <- filter_to_grid(x$filter, tg)
  proxy <- list(time = tg, spikes = x$spikes,
                mu = fg$mu, params = x$params_hat)

  Lk   <- compute_time_rescaling(proxy)
  Lk   <- Lk[is.finite(Lk)]
  uk   <- 1 - exp(-Lk)
  n    <- length(Lk)
  vk   <- qnorm(uk); vk_f <- vk[is.finite(vk)]
  ks_p <- tryCatch(ks.test(uk, "punif")$p.value, error = function(e) NA_real_)
  band <- 1.36 / sqrt(n)

  # ---- (1) uniform QQ ----
  theo <- ppoints(n)
  d1   <- data.frame(theo = theo, emp = sort(uk))
  p1 <- ggplot2::ggplot(d1, ggplot2::aes(x = theo, y = emp)) +
    ggplot2::geom_abline(slope = 1, colour = .COL["warn"], linewidth = 0.75) +
    ggplot2::geom_line(data = data.frame(x = theo, y = pmin(theo + band, 1)),
                       ggplot2::aes(x = x, y = y),
                       colour = .COL["warn"], linetype = "dashed", linewidth = 0.45) +
    ggplot2::geom_line(data = data.frame(x = theo, y = pmax(theo - band, 0)),
                       ggplot2::aes(x = x, y = y),
                       colour = .COL["warn"], linetype = "dashed", linewidth = 0.45) +
    ggplot2::geom_point(size = 0.55, alpha = 0.45, colour = "grey30") +
    ggplot2::labs(x = "Uniform theoretical", y = expression(u[k]),
                  title = "Uniform Q-Q",
                  subtitle = sprintf("KS p = %.3f  (n = %d)", ks_p, n)) +
    .theme_sde_ig()

  # ---- (2) normal-score QQ ----
  nq <- qnorm(ppoints(length(vk_f)))
  d2 <- data.frame(theo = nq, emp = sort(vk_f))
  p2 <- ggplot2::ggplot(d2, ggplot2::aes(x = theo, y = emp)) +
    ggplot2::geom_abline(slope = 1, colour = .COL["warn"], linewidth = 0.75) +
    ggplot2::geom_point(size = 0.55, alpha = 0.45, colour = "grey30") +
    ggplot2::labs(x = "N(0,1) theoretical",
                  y = expression(v[k] == Phi^{-1}(u[k])),
                  title = "Normal-score Q-Q") +
    .theme_sde_ig()

  # ---- (3) ACF ----
  ac_obj <- acf(uk, lag.max = 25, plot = FALSE)
  ci_ac  <- qnorm(0.975) / sqrt(n)
  d3 <- data.frame(lag = as.numeric(ac_obj$lag[-1L]),
                   acf = as.numeric(ac_obj$acf[-1L]))
  p3 <- ggplot2::ggplot(d3, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = c(-ci_ac, ci_ac),
                        linetype = "dashed", colour = .COL["warn"], linewidth = 0.45) +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0),
                          colour = "#0072B2", linewidth = 0.65) +
    ggplot2::geom_point(colour = "#0072B2", size = 1.4) +
    ggplot2::labs(x = "Lag", y = "ACF",
                  title = expression(paste("Residual ACF of ", u[k]))) +
    .theme_sde_ig()

  # ---- (4) Lambda ECDF vs Exp(1) ----
  lk_s <- sort(Lk)
  d4   <- data.frame(x    = lk_s,
                     ecdf = seq_len(n) / (n + 1L),
                     exp1 = 1 - exp(-lk_s))
  p4 <- ggplot2::ggplot(d4) +
    ggplot2::geom_line(ggplot2::aes(x = x, y = exp1),
                       colour = .COL["warn"], linewidth = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = x, y = ecdf),
                        size = 0.55, alpha = 0.4, colour = "grey30") +
    ggplot2::labs(x = expression(Lambda[k]),
                  y = "CDF",
                  title = expression(paste(Lambda[k], " vs Exp(1)")),
                  subtitle = "Orange: theoretical  |  grey: empirical") +
    .theme_sde_ig()

  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2L) +
    patchwork::plot_annotation(
      title    = "Time-Rescaling Goodness-of-Fit Diagnostics",
      subtitle = sprintf("%d inter-beat intervals  |  duration %.1f s",
                         n, tail(x$spikes, 1L)),
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 13),
        plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9)
      )
    )
}
