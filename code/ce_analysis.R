# ---
# Cost-effectiveness analysis
# Author: Darwin Del Castillo
# ---

# Setting up parameters
set.seed(123)
sim_run <- 10000
options(scipen = 999999, max.print = 99999999)

# Function to convert rate to probability
rate_to_prob <- function(r, per, to) {
  1 - exp(-(r / per) * to)
}

# -----------------------------------------------------------
# Parameters
# -----------------------------------------------------------
# Demographic and clinical parameters
age_init <- 42.3  # Mean age from trial
bmi <- 27.2       # Mean BMI from trial # To model the effect of BMI on HTA-HR
sbp_mean <- 113.1 # Mean SBP from trial # To model the effect of BMI on HTA-HR
dbp_mean <- 72    # Mean DBP from trial # To model the effect of BMI on HTA-HR

# Transition probabilities for control group
p_1 <- 0.0546     # Healthy to hypertension (control)

# Intervention case scenario
hr_hta <- 0.49    # Hazard ratio for hypertension with intervention

# Costs
cost_healthy <- 0
cost_hta <- 584
cost_salt <- 1.18
cost_social <- 53.62
cost_intervention <- cost_salt + cost_social

# Cycle correction
cycle_length <- 1  # 1 year

# Utility weights
u_hta <- 0.85      # Mean utility for hypertension (SD ≈ 0.12)

# Discount rate
discount_rate <- 0.03

# Load Peru WHO life expectancy table
# life_expectancy_table <- load_peru_who_life_table()

# -----------------------------------------------------------
# Create transition matrices
# -----------------------------------------------------------
create_transition_matrix <- function(p_1, hr_hta = 1) {
  # Two‑state model: Healthy (state 1) and Hypertension (state 2)
  p_1_adj <- 1 - exp(-(p_1 * hr_hta))   # Adjusted probability Healthy → HTA
  matrix(c(
    # From Healthy
    1 - p_1_adj, p_1_adj,
    # From Hypertension (absorbing)
    0,           1
  ), nrow = 2, byrow = TRUE)
}

# Create matrices
mat_control      <- create_transition_matrix(p_1, hr_hta = 1)
mat_intervention <- create_transition_matrix(p_1, hr_hta = hr_hta)

# -----------------------------------------------------------
# Discounting function
# -----------------------------------------------------------
apply_discount <- function(value, rate, time) {
  value / (1 + rate)^time
}

# -----------------------------------------------------------
# Run deterministic model (two-state version)
# -----------------------------------------------------------
run_deterministic_model <- function(transition_matrix, costs, utility_weights,
                                    cycles, init_state, is_intervention = FALSE,
                                    discount_rate = 0.03) {

  if (is_intervention) {
    base_costs <- c(costs$intervention_healthy, costs$intervention_hta)
  } else {
    base_costs <- c(costs$healthy, costs$hta)
  }

  util_vector <- c(1, utility_weights$hta)

  state_dist <- matrix(0, nrow = cycles + 1, ncol = 2)
  state_dist[1, ] <- init_state

  costs_by_cycle <- numeric(cycles)
  qaly_by_cycle   <- numeric(cycles)

  for (t in 1:cycles) {
    if (t > 1) state_dist[t, ] <- state_dist[t - 1, ] %*% transition_matrix

    cycle_costs <- sum(state_dist[t, ] * base_costs)
    cycle_qaly  <- sum(state_dist[t, ] * util_vector)

    costs_by_cycle[t] <- cycle_costs / (1 + discount_rate)^(t - 1)
    qaly_by_cycle[t]  <- cycle_qaly / (1 + discount_rate)^(t - 1)
  }

  total_costs <- sum(costs_by_cycle)
  total_qalys <- sum(qaly_by_cycle)

  list(
    state_dist     = state_dist,
    costs_by_cycle = costs_by_cycle,
    qaly_by_cycle  = qaly_by_cycle,
    total_costs    = total_costs,
    total_qalys    = total_qalys
  )
}

# -----------------------------------------------------------
# Run both strategies
# -----------------------------------------------------------
costs <- list(
  healthy      = cost_healthy,
  hta          = cost_hta,
  intervention_healthy = cost_healthy + cost_intervention,
  intervention_hta = cost_hta + cost_intervention
)

utility_weights <- list(
  hta = u_hta
)

init_states <- c(1, 0)  # All start in Healthy state
cycles <- 60  # Lifetime horizon

# Run control strategy
control_results <- run_deterministic_model(
  transition_matrix = mat_control,
  costs = costs,
  utility_weights = utility_weights,
  cycles = cycles,
  init_state = init_states,
  is_intervention = FALSE
)

# Run intervention strategy
intervention_results <- run_deterministic_model(
  transition_matrix = mat_intervention,
  costs = costs,
  utility_weights = utility_weights,
  cycles = cycles,
  init_state = init_states,
  is_intervention = TRUE
)

# -----------------------------------------------------------
# Calculate incremental cost-effectiveness ratio (ICER)
# -----------------------------------------------------------
incremental_cost <- intervention_results$total_costs - control_results$total_costs
incremental_effect <- intervention_results$total_qalys - control_results$total_qalys
icer <- incremental_cost / incremental_effect      # cost per QALY gained

# -----------------------------------------------------------
# Summarize results
# -----------------------------------------------------------
summary_results <- data.frame(
  Strategy         = c("Control", "Intervention"),
  Total_Costs      = c(control_results$total_costs, intervention_results$total_costs),
  Total_QALYs      = c(control_results$total_qalys, intervention_results$total_qalys),
  Incremental_Cost = c(NA, incremental_cost),
  Incremental_Effect = c(NA, incremental_effect),
  ICER             = c(NA, icer)
)

print(summary_results)

# -----------------------------------------------------------
# Probabilistic Sensitivity Analysis (PSA)
# -----------------------------------------------------------
library(parallel)
library(ggplot2)

# Number of PSA simulations
N <- sim_run

# Function to draw one PSA sample and run the deterministic model
run_psa_iteration <- function(i) {
  # Sample demographic and clinical parameters
  age_init_i <- rnorm(1, mean = 43.3, sd = 17.2 / sqrt(2376)) # nolint: object_usage_linter.
  bmi_i <- rnorm(1, mean = 27.2, sd = 4.6 / sqrt(2376)) # nolint: object_usage_linter.
  sbp_sd <- 17.0 / sqrt(2376)
  sbp_meanlog <- log(113.1) - 0.5 * (sbp_sd / 113.1)^2
  sbp_sdlog <- sqrt(log(1 + (sbp_sd / 113.1)^2))
  sbp_mean_i <- rlnorm(1, meanlog = sbp_meanlog, sdlog = sbp_sdlog) # nolint: object_usage_linter.
  dbp_sd <- 10.1 / sqrt(2376)
  dbp_meanlog <- log(72) - 0.5 * (dbp_sd / 72)^2
  dbp_sdlog <- sqrt(log(1 + (dbp_sd / 72)^2))
  dbp_mean_i <- rlnorm(1, meanlog = dbp_meanlog, sdlog = dbp_sdlog) # nolint: object_usage_linter.

  # Hazard ratio for hypertension
  hr_sdlog <- (log(0.66) - log(0.31)) / 3.92
  hr_hta_i <- rlnorm(1, meanlog = log(0.45), sdlog = hr_sdlog)

  # Costs
  cost_hta_i <- rgamma(1, shape = 476, rate = 1)
  # Fixed intervention costs remain the same as defined globally
  costs_i <- list(
    healthy      = cost_healthy,
    hta          = cost_hta_i,
    intervention_healthy = cost_healthy + cost_intervention,
    intervention_hta = cost_hta_i + cost_intervention
  )
  # Utility weights
  # Beta(6.67, 1.18) gives mean 0.85 and sd ≈ 0.12
  u_hta_i <- rbeta(1, shape1 = 6.67, shape2 = 1.18)
  utils_i <- list(hta = u_hta_i)

  # Create transition matrices with sampled hr
  mat_control_i      <- create_transition_matrix(p_1, hr_hta = 1)
  mat_intervention_i <- create_transition_matrix(p_1, hr_hta = hr_hta_i)

  # Run deterministic model for both strategies
  res_control_i <- run_deterministic_model(
    transition_matrix = mat_control_i,
    costs = costs_i,
    utility_weights = utils_i,
    cycles = cycles,
    init_state = init_states,
    is_intervention = FALSE
  )
  res_int_i <- run_deterministic_model(
    transition_matrix = mat_intervention_i,
    costs = costs_i,
    utility_weights = utils_i,
    cycles = cycles,
    init_state = init_states,
    is_intervention = TRUE
  )

  # Compute incremental cost and effect
  inc_cost <- res_int_i$total_costs - res_control_i$total_costs
  inc_eff <- res_int_i$total_qalys - res_control_i$total_qalys

  return(c(inc_cost = inc_cost, inc_eff = inc_eff))
}

# Run PSA in parallel
nCores <- detectCores() - 1
psa_out <- mclapply(1:N, run_psa_iteration, mc.cores = nCores)
psa_df <- do.call(rbind, psa_out)
psa_df <- as.data.frame(psa_df)
psa_df$icer <- psa_df$inc_cost / psa_df$inc_eff

# Summarize PSA ICER
icer_mean <- mean(psa_df$icer, na.rm = TRUE)
icer_ci <- quantile(psa_df$icer, c(0.025, 0.975), na.rm = TRUE)
cat("PSA ICER mean:", round(icer_mean, 2), "95% CI [", round(icer_ci[1], 2), ",", round(icer_ci[2], 2), "]\n")

# -----------------------------------------------------------
# Graphs for PSA
# -----------------------------------------------------------

# Peruvian QALY threshold: S/25,007 per AVAC converted to 2024 USD
# Average 2024 USD/PEN exchange rate = 3.7538
wtp_threshold <- 25007 / 3.7538  # ≈ 6662 USD per QALY
# Cost-effectiveness plane
ce_plane_df <- psa_df
# Calculate centroid
centroid_df <- data.frame(
  inc_eff = mean(ce_plane_df$inc_eff, na.rm = TRUE),
  inc_cost = mean(ce_plane_df$inc_cost, na.rm = TRUE)
)
gg_ce_plane <- ggplot(ce_plane_df, aes(x = inc_eff, y = inc_cost)) +
  geom_point(color = "#b20000a3", alpha = 0.4) +
  geom_abline(intercept = 0, slope = wtp_threshold, color = "#0d009e", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 0.6) +
  geom_hline(yintercept = 0, color = "black", size = 0.6) +
  stat_ellipse(type = "norm", level = 0.95, color = "#999999", linetype = "solid", size = 1) +
  geom_point(data = centroid_df, aes(x = inc_eff, y = inc_cost), color = "#999999", size = 4, shape = 8) +
  labs(
    x = "Incremental QALYs Gained",
    y = "Incremental Costs",
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
print(gg_ce_plane)

ggsave(gg_ce_plane,
filename = "Cost effectiveness plane.jpeg",
width = 20, height = 40, units = "cm",
dpi = 800,
path = "output/figs")

# Cost-effectiveness acceptability curve (CEAC)
wtp_vals <- seq(0, 10000, by = 100)
ceac <- sapply(wtp_vals, function(wtp) {
  mean(psa_df$inc_cost <= wtp * psa_df$inc_eff, na.rm = TRUE)
})
ceac_df <- data.frame(
  wtp = wtp_vals,
  probability = ceac
)
gg_ceac <- ggplot(ceac_df, aes(x = wtp, y = probability)) +
  geom_line(color = "#bf0000", size = 1) +
  labs(
    x = "Willingness-to-pay Threshold",
    y = "Probability Cost-effective",
    title = "Cost-effectiveness Acceptability Curve"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
print(gg_ceac)

# Expected Value of Perfect Information (EVPI)
evpi <- sapply(wtp_vals, function(wtp) {
  # Net monetary benefit for the intervention
  nb_int <- wtp * psa_df$inc_eff - psa_df$inc_cost
  # Control NB is always zero

  # 1) Expected value with perfect information: E[max(NB_i)]
  ev_with_info <- mean(pmax(nb_int, 0))

  # 2) Expected value under current decision: max(E[NB_i])
  ev_current <- max(mean(nb_int), 0)

  # 3) EVPI = E[max(NB)] - max(E[NB])
  ev_with_info - ev_current
})

# Build data.frame and plot
evpi_df <- data.frame(
  wtp  = wtp_vals,
  evpi = evpi
)

gg_evpi <- ggplot(evpi_df, aes(x = wtp, y = evpi)) +
  geom_line() +
  labs(
    x     = "Willingness-to-pay Threshold",
    y     = "EVPI",
    title = "Expected Value of Perfect Information"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(gg_evpi)