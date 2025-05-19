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

# ------------------------------------------------------------
# Function to load the Peru WHO life expectancy table
# ------------------------------------------------------------
# Life expectancy from WHO for Peru in 2019
load_peru_who_life_table <- function() {
# Life expectancy table for Peru
  peru_life_table <- data.frame(
    age_group = c("<1 year", "1-4 years", "5-9 years", "10-14 years", "15-19 years", 
                  "20-24 years", "25-29 years", "30-34 years", "35-39 years",
                  "40-44 years", "45-49 years", "50-54 years", "55-59 years",
                  "60-64 years", "65-69 years", "70-74 years", "75-79 years",
                  "80-84 years", "85+ years"),
    both_sexes = c(79.8978, 79.6661, 75.9019, 71.0033, 66.1034, 61.2717, 56.503, 
                   51.7672, 47.0478, 42.3573, 37.7193, 33.1538, 28.6979, 24.3816, 
                   20.2588, 16.3779, 12.7774, 9.5738, 6.6763),
    male = c(78.4557, 78.276, 74.5268, 69.6348, 64.7382, 59.9354, 55.2476,
             50.6094, 45.98, 41.3695, 36.8012, 32.2972, 27.906, 23.6606,
             19.6184, 15.8347, 12.3418, 9.2986, 6.5418),
    female = c(81.3357, 81.0453, 77.2645, 72.359, 67.4556, 62.5952, 57.7503,
               52.9122, 48.0938, 43.3145, 38.6001, 33.9677, 29.4436, 25.0526,
               20.8449, 16.8634, 13.1532, 9.7974, 6.7664)
  )
  
  # Create vector to store life expectancy by individual age
  max_age <- 100 #Maximum age in the trial
  life_expectancy <- numeric(max_age + 1)
  
  # Function to convert age group to numeric range
  get_age_range <- function(age_group) {
    if (age_group == "<1 year") {
      return(list(min = 0, max = 0))
    } else if (age_group == "85+ years") {
      return(list(min = 85, max = max_age))
    } else {
      # Extract numbers from age group (e.g. "1-4 years" -> min=1, max=4)
      range_str <- gsub(" years", "", age_group)
      range_parts <- strsplit(range_str, "-")[[1]]
      min_age <- as.numeric(range_parts[1])
      max_age <- as.numeric(range_parts[2])
      return(list(min = min_age, max = max_age))
    }
  }
  
  # Fill the life expectancy vector for each individual year of age
  for (i in seq_len(nrow(peru_life_table))) {
    age_range <- get_age_range(peru_life_table$age_group[i])
    for (age in age_range$min:age_range$max) {
      life_expectancy[age + 1] <- peru_life_table$both_sexes[i]
    }
  }
  
  return(life_expectancy)
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
hr_hta <- 0.45    # Hazard ratio for hypertension with intervention

# Costs
cost_healthy <- 0
cost_hta <- 584
cost_intervention <- 1.18

# Cycle correction
cycle_length <- 1  # 1 year

# Disability weights
d_hta <- 0.02      # DW for hypertension

# Discount rate
discount_rate <- 0.03

# Load Peru WHO life expectancy table
life_expectancy_table <- load_peru_who_life_table()

# -----------------------------------------------------------
# Create transition matrices
# -----------------------------------------------------------
create_transition_matrix <- function(p_1, hr_hta = 1) {
  # Two‑state model: Healthy (state 1) and Hypertension (state 2)
  p_1_adj <- 1 - exp(-(p_1 * hr_hta))        # Adjusted probability Healthy → HTA
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
run_deterministic_model <- function(transition_matrix, costs, disability_weights,
                                    cycles, init_state, is_intervention = FALSE,
                                    discount_rate = 0.03) {

  base_costs <- c(costs$healthy, costs$hta)
  if (is_intervention) {
    base_costs <- base_costs + c(costs$intervention, costs$intervention)
  }

  dw_vector <- c(0, disability_weights$hta)

  state_dist <- matrix(0, nrow = cycles + 1, ncol = 2)
  state_dist[1, ] <- init_state

  costs_by_cycle <- numeric(cycles)
  yld_by_cycle   <- numeric(cycles)

  for (t in 1:cycles) {
    if (t > 1) state_dist[t, ] <- state_dist[t - 1, ] %*% transition_matrix

    cycle_costs <- sum(state_dist[t, ] * base_costs)
    cycle_yld   <- sum(state_dist[t, ] * dw_vector)

    costs_by_cycle[t] <- cycle_costs / (1 + discount_rate)^(t - 1)
    yld_by_cycle[t]   <- cycle_yld   / (1 + discount_rate)^(t - 1)
  }

  total_costs <- sum(costs_by_cycle)
  total_dalys <- sum(yld_by_cycle)  # DALYs = YLD only in two‑state model

  list(
    state_dist     = state_dist,
    costs_by_cycle = costs_by_cycle,
    yld_by_cycle   = yld_by_cycle,
    total_costs    = total_costs,
    total_dalys    = total_dalys
  )
}

# -----------------------------------------------------------
# Run both strategies
# -----------------------------------------------------------
costs <- list(
  healthy      = cost_healthy,
  hta          = cost_hta,
  intervention = cost_intervention
)

disability_weights <- list(
  hta = d_hta
)

init_states <- c(1, 0)  # All start in Healthy state
cycles <- 60  # Lifetime horizon

# Run control strategy
control_results <- run_deterministic_model(
  transition_matrix = mat_control,
  costs = costs,
  disability_weights = disability_weights,
  cycles = cycles,
  init_state = init_states,
  is_intervention = FALSE
)

# Run intervention strategy
intervention_results <- run_deterministic_model(
  transition_matrix = mat_intervention,
  costs = costs,
  disability_weights = disability_weights,
  cycles = cycles,
  init_state = init_states,
  is_intervention = TRUE
)

# -----------------------------------------------------------
# Calculate incremental cost-effectiveness ratio (ICER)
# -----------------------------------------------------------
incremental_cost <- intervention_results$total_costs - control_results$total_costs
incremental_effect <- control_results$total_dalys - intervention_results$total_dalys
icer <- incremental_cost / incremental_effect

# -----------------------------------------------------------
# Summarize results
# -----------------------------------------------------------
summary_results <- data.frame(
  Strategy           = c("Control", "Intervention"),
  Total_Costs        = c(control_results$total_costs, intervention_results$total_costs),
  Total_DALYs        = c(control_results$total_dalys, intervention_results$total_dalys),
  Incremental_Cost   = c(NA, incremental_cost),
  Incremental_Effect = c(NA, incremental_effect),
  ICER               = c(NA, icer)
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
  age_init_i <- rnorm(1, mean = 43.3, sd = 17.2 / sqrt(2376))
  bmi_i <- rnorm(1, mean = 27.2, sd = 4.6 / sqrt(2376))
  sbp_sd <- 17.0 / sqrt(2376)
  sbp_meanlog <- log(113.1) - 0.5 * (sbp_sd / 113.1)^2
  sbp_sdlog <- sqrt(log(1 + (sbp_sd / 113.1)^2))
  sbp_mean_i <- rlnorm(1, meanlog = sbp_meanlog, sdlog = sbp_sdlog)
  dbp_sd <- 10.1 / sqrt(2376)
  dbp_meanlog <- log(72) - 0.5 * (dbp_sd / 72)^2
  dbp_sdlog <- sqrt(log(1 + (dbp_sd / 72)^2))
  dbp_mean_i <- rlnorm(1, meanlog = dbp_meanlog, sdlog = dbp_sdlog)

  # Hazard ratio for hypertension
  hr_sdlog <- (log(0.66) - log(0.31)) / 3.92
  hr_hta_i <- rlnorm(1, meanlog = log(0.45), sdlog = hr_sdlog)

  # Costs
  cost_hta_i <- rgamma(1, shape = 476, rate = 1)
  # Fixed intervention costs remain the same as defined globally
  costs_i <- list(
    healthy      = cost_healthy,
    hta          = cost_hta_i,
    intervention = cost_intervention
  )
  # Disability weights
  d_hta_i <- rbeta(1, 2, 98)
  dws_i <- list(hta = d_hta_i)

  # Create transition matrices with sampled hr
  mat_control_i      <- create_transition_matrix(p_1, hr_hta = 1)
  mat_intervention_i <- create_transition_matrix(p_1, hr_hta = hr_hta_i)

  # Run deterministic model for both strategies
  res_control_i <- run_deterministic_model(
    transition_matrix = mat_control_i,
    costs = costs_i,
    disability_weights = dws_i,
    cycles = cycles,
    init_state = init_states,
    is_intervention = FALSE
  )
  res_int_i <- run_deterministic_model(
    transition_matrix = mat_intervention_i,
    costs = costs_i,
    disability_weights = dws_i,
    cycles = cycles,
    init_state = init_states,
    is_intervention = TRUE
  )

  # Compute incremental cost and effect
  inc_cost <- res_int_i$total_costs - res_control_i$total_costs
  inc_eff <- res_control_i$total_dalys - res_int_i$total_dalys

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

# Cost-effectiveness plane
ce_plane_df <- psa_df
# Calculate centroid
centroid_df <- data.frame(
  inc_eff = mean(ce_plane_df$inc_eff, na.rm = TRUE),
  inc_cost = mean(ce_plane_df$inc_cost, na.rm = TRUE)
)
gg_ce_plane <- ggplot(ce_plane_df, aes(x = inc_eff, y = inc_cost)) +
  geom_point(color = "#0072B2", alpha = 0.4) +
  stat_ellipse(type = "norm", level = 0.95, color = "#999999", linetype = "dashed", size = 1) +
  geom_point(data = centroid_df, aes(x = inc_eff, y = inc_cost), color = "#D55E00", size = 4, shape = 8) +
  labs(
    x = "Incremental DALYs Averted",
    y = "Incremental Costs",
    title = "Cost-effectiveness Plane"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
print(gg_ce_plane)

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
  theme_minimal()

print(gg_evpi)