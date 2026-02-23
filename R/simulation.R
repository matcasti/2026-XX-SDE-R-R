# ==============================================================================
# Inverse Gaussian Point Process Model for HRV
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SIMULATION ENGINE
# ------------------------------------------------------------------------------

#' Simulate Heartbeat Series with SDE-Driven Autonomic Dynamics
#'
#' @param duration_secs Total time to simulate (seconds)
#' @param dt Simulation time step (high resolution recommended, e.g., 0.005)
#' @param stress_fn A function taking time 't' and returning scalar input u(t)
#' @return List containing time grid, state trajectories, and spike times
simulate_heart_sde_ig <- function(duration_secs, dt = 0.005, stress_fn = NULL) {

  # Default stressor: No stress if function not provided
  if(is.null(stress_fn)) stress_fn <- function(t) { 0 }

  times <- seq(0, duration_secs, by = dt)
  n <- length(times)

  # --- Model Parameters (The "Physiology") ---
  # Parasympathetic (Fast decay, High noise)
  a_p <- 0.4; sig_p <- 0.2; k_par <- 0.2
  # Sympathetic (Slow decay, Low noise)
  a_s <- 0.04; sig_s <- 0.1; k_sym <- 0.4
  # Coupling effects
  a_sp <- -0.05; a_ps <- -0.1;

  # Pacemaker properties
  rho_0 <- 1.1; kappa <- 20

  # --- Initialization ---
  s <- numeric(n); p <- numeric(n)
  u <- numeric(n)
  spikes <- numeric(0); last_spike_t <- -2.0

  # --- Main Loop ---
  for(i in 1:(n-1)) {
    t_curr <- times[i]
    u_val <- stress_fn(t_curr)
    u[i] <- u_val

    # 1. Euler-Maruyama Update for SDEs
    # p(t) is suppressed by stress, s(t) is excited by stress
    dp <- (-a_p * p[i] + a_ps * s[i] - 0.4 * u_val) * dt + sig_p * rnorm(1, 0, sqrt(dt))
    ds <- (-a_s * s[i] + a_sp * p[i] + 0.1 * u_val) * dt + sig_s * rnorm(1, 0, sqrt(dt))

    p[i+1] <- p[i] + dp
    s[i+1] <- s[i] + ds

    # 2. Map States to Mean Interval (mu)
    current_rate <- rho_0 + k_sym * s[i+1] - k_par * p[i+1]
    current_rate <- max(0.3, current_rate) # Floor at ~24 BPM
    mu <- 1 / current_rate

    # 3. Calculate Inverse Gaussian Hazard

    ## Time since last heartbeat
    tau <- t_curr - last_spike_t

    if(tau < 0.02) {
      lambda_val <- 0 # Refractory period
    } else {
      # Numerical safeguards for large exponents
      exp_arg <- 2 * kappa / mu
      term_cdf_part <- exp(exp_arg) * pnorm(-sqrt(kappa/tau)*(tau/mu + 1))
      if(exp_arg > 700) {
        # If huge, CDF term explodes/saturates
        term_cdf_part <- 0
      } else {
        term_cdf_part <- exp(exp_arg) * pnorm(-sqrt(kappa/tau)*(tau/mu + 1))
      }

      term_pdf <- sqrt(kappa/(2*pi*tau^3)) * exp(-kappa*(tau-mu)^2/(2*mu^2*tau))
      term_cdf <- pnorm(sqrt(kappa/tau)*(tau/mu - 1)) + term_cdf_part

      # Hazard = PDF / Survival
      lambda_val <- term_pdf / max(1e-9, (1 - term_cdf))
      if(is.na(lambda_val) || is.infinite(lambda_val)) lambda_val <- 0
    }

    # 4. Bernoulli Trial (Thinning)
    if(runif(1) < lambda_val * dt) {
      spikes <- c(spikes, times[i+1])
      last_spike_t <- times[i+1]
    }
  }

  return(list(time=times, s=s, p=p, spikes=spikes, u=u))
}

# ------------------------------------------------------------------------------
# 2. MAIN EXECUTION BLOCK
# ------------------------------------------------------------------------------

# --- A. Generate Synthetic Data ---
cat("\n--- 1. Generating Synthetic 'Patient' Data ---\n")

stress_scenario <- function(t) {
  if(t < 300) return(0)
  if(t >= 300 & t < 360) return((t-300)/60) # Linear ramp 0->1
  if(t >= 360 & t < 420) return(1) # 1
  return(0) # Recovery
}

set.seed(123)
sim_data <- simulate_heart_sde_ig(duration_secs = 720, stress_fn = stress_scenario)

# Extract just the "observed" spike times for the filter
observed_spikes <- sim_data$spikes
cat(paste("Generated", length(observed_spikes), "heartbeats.\n"))

# --- C. Visualization ---
cat("\n--- 3. Visualizing Results ---\n")

# Layout: 3 Plots stackedX
par(mfrow=c(3,1), mar=c(3,4,2,2), oma=c(2,0,0,0))

# Plot 1: R-R Intervals (Tachogram)
rr_intervals <- diff(observed_spikes)
rr_times <- observed_spikes[-1]
plot(rr_times, rr_intervals, type='l', pch=19, cex=0.6, col="grey30",
     ylab="R-R Interval (s)", main="Observed Heart Rate Variability", frame.plot=F,
     ylim = c(0, 2))
grid()

# Plot 2: Parasympathetic Recovery
plot(sim_data$time, sim_data$p, type='l', col="#0072B2", lwd=1,
     ylab="Parasympathetic Tone", main="Vagal Tone Estimation",
     frame.plot=F)
grid()

# Plot 3: Sympathetic Recovery
plot(sim_data$time, sim_data$s, type='l', col="#D55E00", lwd=1,
     ylab="Sympathetic Tone", main="Sympathetic Tone Estimation",
     frame.plot=F)
grid()
