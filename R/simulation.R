# sim_and_viz.R
# Simulate 2-latent SDEs (PNS = p, SNS = s) and Poisson-bin heartbeat counts.
# Visualize latents, smoothed HR (bpm), and beat raster.

set.seed(123)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)

# ---- simulation settings ----
total_minutes <- 12
T_sec <- total_minutes * 60   # 720 s
dt <- 0.1                    # time step (s)
time <- seq(0, T_sec, by = dt)
K <- length(time)

# protocol: 5 min rest, 2 min exercise, 5 min rest
u <- rep(0, K)
u[time >= 5*60 & time < 7*60] <- 1
# optional: extend easily for 24h by changing total_minutes

# ---- true parameters (physiological-ish) ----
a_p <- 1/2        # PNS decay rate => tau_p = 2 s
a_s <- 1/60       # SNS decay rate => tau_s = 60 s
b_ps <- -0.02     # SNS -> PNS coupling (suppression)
b_sp <- -0.02    # PNS -> SNS coupling (small)
c_p <- -0.8       # input gain on PNS (exercise reduces PNS)
c_s <-  1.0       # input gain on SNS (exercise increases SNS)
sigma_p <- 0.25
sigma_s <- 0.25

# observation params
mu <- 0.0         # baseline log-rate (log beats per second) => exp(0)=1 Hz = 60 bpm
alpha_p <- 1.0
alpha_s <- 1.0
kappa <- -3.0     # simple 1-lag history effect (strong refractoriness)
# note: history term multiplies previous bin count (0/1 typical for small dt)

# ---- simulate latent SDEs (Euler-Maruyama) ----
p <- numeric(K); s <- numeric(K)
p[1] <- 0; s[1] <- 0  # initial states
for (k in 1:(K-1)) {
  dp_det <- (-a_p * p[k] + b_ps * s[k] + c_p * u[k]) * dt
  ds_det <- (-a_s * s[k] + b_sp * p[k] + c_s * u[k]) * dt
  p[k+1] <- p[k] + dp_det + rnorm(1, 0, sigma_p * sqrt(dt))
  s[k+1] <- s[k] + ds_det + rnorm(1, 0, sigma_s * sqrt(dt))
}

# ---- generate counts by Poisson approximation over each dt bin ----
lambda <- exp(mu + alpha_p * p - alpha_s * s)  # instantaneous rate in beats/sec
# incorporate 1-lag history: implement sequentially
y <- integer(K)
for (k in 1:K) {
  hist_term <- ifelse(k==1, 0, kappa * y[k-1])
  rate_k <- exp(mu + alpha_p * p[k] - alpha_s * s[k] + hist_term)
  # Poisson with mean rate * dt
  y[k] <- rpois(1, lambda = rate_k * dt)
  # avoid more than 1 per tiny dt for visualization realism (optional)
  # if (y[k] > 1) y[k] <- 1
}

# reconstruct beat times (place counts at bin centers, if count>0 put that many ticks)
beat_times <- rep(time, times = y)

# ---- compute smoothed HR (bpm) for plotting: moving window counts per 60s -> bpm ----
# instantaneous rate estimate via moving sum over 5-second window, scaled to bpm
window_s <- 5       # seconds
window_bins <- round(window_s / dt)
counts_ts <- y
smoothed_rate_bps <- rollapply(counts_ts, width = window_bins, FUN = function(x) sum(x) / (window_s), align="right", fill=NA)
smoothed_bpm <- smoothed_rate_bps * 60

# package for plotting
df <- tibble(
  t = time,
  p = p,
  s = s,
  lambda = lambda,
  count = y,
  smoothed_bpm = smoothed_bpm,
  u = u
)

plot(smoothed_bpm ~ t, type = "l", df)

# ---- visualization ----
# 1) Latent states + protocol
p1 <- ggplot(df, aes(x = t)) +
  annotate(geom = "rect", xmin = 5*60, xmax = 7*60, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.3) +
  geom_line(aes(y = p, color = "PNS (p)")) +
  geom_line(aes(y = s, color = "SNS (s)")) +
  scale_color_manual("", values = c("PNS (p)" = "steelblue", "SNS (s)" = "firebrick")) +
  labs(x = "Time (s)", y = "Latent state (arb. units)", title = "Simulated latent states (PNS & SNS)\norange = exercise window") +
  theme_minimal()

# 2) Smoothed HR (bpm) + true instantaneous rate (lambda*60)
p2 <- ggplot(df, aes(x = t)) +
  geom_line(aes(y = smoothed_bpm), linewidth = 0.9, na.rm = TRUE) +
  geom_line(aes(y = lambda * 60), linetype = "dashed", alpha = 0.7) +
  labs(x = "Time (s)", y = "Heart rate (bpm)", title = "Smoothed HR (window 5s) and true instantaneous HR (dashed)") +
  theme_minimal()

# 3) raster of beat times + counts
if (length(beat_times) > 0) {
  beat_df <- tibble(t = beat_times, id = 1:length(beat_times))
  p3 <- ggplot() +
    geom_point(data = beat_df, aes(x = t, y = id %% 50), shape = '|', size = 3) +
    labs(x = "Time (s)", y = "Beat raster (mod 50)", title = "Simulated beat raster (tick per beat)") +
    theme_minimal()
} else {
  p3 <- ggplot() + ggtitle("No beats simulated (unexpected).")
}

# show plots
print(p1)
print(p2)
print(p3)

# save plots optionally
# ggsave("latents.png", p1, width=10, height=4)
# ggsave("hr.png", p2, width=10, height=3)
# ggsave("raster.png", p3, width=10, height=2.5)
