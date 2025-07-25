
# ---
# Missing Parameters Estimation
# Author: Darwin Del Castillo
# ---
# ----------------------------------------------------------
# Code for estimating missing parameters for the MCMC model
# ----------------------------------------------------------

## CVD events due to pre-HTA
### Inputs
events_no_htn <- 8               # from Lazo‑Porras et al. (2022)
person_time_no_htn <- 12789      # person‑years
rr_pre_mean <- 1.55              # pooled RR from Huang et al. (2013)
rr_pre_low  <- 1.41              # 95 % CI lower
rr_pre_high <- 1.71              # 95 % CI upper
n_iter <- 10000                  # Monte‑Carlo draws

### Draw lambda_0 (normotensive incidence)
shape_gamma <- events_no_htn
rate_gamma <- person_time_no_htn             # gamma(rate=…) uses 1/scale
lambda0_draw <- rgamma(n_iter,
                       shape = shape_gamma,
                       rate  = rate_gamma)   # per person‑year

### Draw RR_pre (log‑normal)
sigma_logn <- (log(rr_pre_high) - log(rr_pre_low)) / (2 * 1.96)
mu_logn    <- log(rr_pre_mean)

rr_pre_draw <- rlnorm(n_iter,
                      meanlog = mu_logn,
                      sdlog   = sigma_logn)

### Derived incidence for pre‑HTA
lambda_pre_draw <- lambda0_draw * rr_pre_draw   # per person‑year

### Summaries
ci <- quantile(lambda_pre_draw, probs = c(0.025, 0.5, 0.975))
ci_100k <- ci * 1e5

print(round(ci, 6))         # per person‑year
print(round(ci_100k, 1))    # per 100 000 person‑years

# ----------------------------------------------------------
# Composite CV event disability weight (mean of stroke and AMI)
# ----------------------------------------------------------

library(truncnorm)

## GBD-2021 DW inputs

dw <- data.frame(
  condition = c("Stroke_mild", "Stroke_mod", "Stroke_Severe", "AMI"),
  mean  = c(0.019, 0.07, 0.552, 0.432),
  low95 = c(0.01, 0.046, 0.376, 0.288),
  hi95  = c(0.032, 0.099, 0.706, 0.579)
)

### derive SD from 95 % UI  (approx - assumes ±1.96 SD either side)
dw$sd <- (dw$hi95 - dw$low95) / (2 * 1.96)

### Monte-Carlo propagation
set.seed(965006311)
n_draw <- 1e5
draws <- sapply(1:nrow(dw), function(i) { # nolint
  rtruncnorm(n_draw, a = 0, b = 1,
             mean = dw$mean[i],
             sd   = dw$sd[i])
})
colnames(draws) <- dw$condition

### mean & variance of composite draws
mu_comp  <- mean(draws)
var_comp <- var(draws)

### 95 % UI
ui_comp  <- quantile(draws, c(0.025, 0.975))


### Fit Beta(α,β) by method of moments
alpha <- mu_comp * ((mu_comp * (1 - mu_comp) / var_comp) - 1)
beta  <- (1 - mu_comp) * ((mu_comp * (1 - mu_comp) / var_comp) - 1)

### Output

list(
  DW_mean = mu_comp,
  DW_95UI = ui_comp,
  Beta_alpha = alpha,
  Beta_beta  = beta
)

# --------------------------------------------------------------------------
# Beta distribution parameters for average disability weight
# --------------------------------------------------------------------------

# Calculate average DW and its approximate CI from the GBD inputs
mean_dw   <- mean(dw$mean)
low95_dw  <- mean(dw$low95)
hi95_dw   <- mean(dw$hi95)

# Derive standard deviation from the 95% UI (assumes ±1.96 SD)
sd_dw     <- (hi95_dw - low95_dw) / (2 * 1.96)
var_dw    <- sd_dw^2

# Fit Beta(α, β) by method of moments for the average DW
alpha_dw  <- mean_dw * ((mean_dw * (1 - mean_dw) / var_dw) - 1)
beta_dw   <- (1 - mean_dw) * ((mean_dw * (1 - mean_dw) / var_dw) - 1)

# Print or return the parameters
list(
  Avg_DW_mean  = mean_dw,
  Avg_DW_SD    = sd_dw,
  Beta_alpha   = alpha_dw,
  Beta_beta    = beta_dw
)

# --------------------------------------------------------------------------
# Pre-hypertension disability weight and hypertension among prehypertensives
# --------------------------------------------------------------------------
## Calculating person-time at risk
### Order data by id and visit
data <- data[order(data$id, data$visit), ]

### Calculate time between visits
data$time_diff <- c(NA, diff(as.numeric(data$visit_date)))
data$time_diff[!duplicated(data$id)] <- as.numeric(data$visit_date[!duplicated(data$id)] - 
                                                    data$date[!duplicated(data$id)])

### Calculate cumulative person-time
data$cum_time <- ave(data$time_diff, data$id, FUN = cumsum)

# Calculate start and stop times for each interval
data$tstart <- c(0, data$cum_time[-nrow(data)])
data$tstop <- data$cum_time

return(data)