
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
# Composite CV event disability weight (moderate severity)
# ----------------------------------------------------------

library(truncnorm)

## GBD-2021 DW inputs (moderate, chronic phase)

dw <- data.frame(
  condition = c("Stroke_mod", "IHD_angina_mod", "HF_mod"),
  mean  = c(0.0701, 0.0795, 0.0717),
  low95 = c(0.0464, 0.0519, 0.0466),
  hi95  = c(0.0993, 0.1129, 0.1031)
)

### derive SD from 95 % UI  (approx - assumes ±1.96 SD either side)
dw$sd <- (dw$hi95 - dw$low95) / (2 * 1.96)

### Peruvian weights
w_raw <- c(Stroke_mod = 0.22,
           IHD_angina_mod = 0.38,
           HF_mod = 0.25)
w <- w_raw / sum(w_raw)   # rescale to sum = 1

### Monte-Carlo propagation
set.seed(123)
n_draw <- 1e5
draws <- sapply(1:nrow(dw), function(i) { # nolint
  rtruncnorm(n_draw, a = 0, b = 1,
             mean = dw$mean[i],
             sd   = dw$sd[i])
})
colnames(draws) <- dw$condition
cv_comp_draw <- draws %*% w

### mean & variance of composite draws
mu_comp  <- mean(cv_comp_draw)
var_comp <- var(cv_comp_draw)

### 95 % UI
ui_comp  <- quantile(cv_comp_draw, c(0.025, 0.975))


### Fit Beta(α,β) by method of moments
alpha <- mu_comp * ((mu_comp * (1 - mu_comp) / var_comp) - 1)
beta  <- (1 - mu_comp) * ((mu_comp * (1 - mu_comp) / var_comp) - 1)

### Output
list(
  Composite_DW_mean = mu_comp,
  Composite_DW_95UI = ui_comp,
  Beta_alpha = alpha,
  Beta_beta  = beta
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