# ---
# Cost-effectiveness analysis
# Author: Darwin Del Castillo
# ---

# Setting up parameters
set.seed(965006311)
sim_run <- 10000
options(scipen = 999999, max.print = 99999999)

# Functions
# Include:
# - Functions to transform life expectancy
# - Functions to transform rates to probabilities

source("functions/functions.R")

# -----------------------------------------------------------
# Parameters
# -----------------------------------------------------------
# Demographic and clinical parameters
age_init <- 43    # Mean age from trial
bmi <- 27.2       # Mean BMI from trial # To model the effect of BMI on HTA-HR
sbp_mean <- 113.1 # Mean SBP from trial # To model the effect of BMI on HTA-HR
dbp_mean <- 72    # Mean DBP from trial # To model the effect of BMI on HTA-HR

# Logistic penalization function for BMI
L_bmi   <- 0.1   # Minimum fraction of effect retained at high BMI
bmi_mid <- 30    # BMI at which half the effect is retained
s_bmi   <- 5     # Slope parameter for transition
frac_bmi_retained <- function(imc) {
  L_bmi + (1 - L_bmi) / (1 + exp((imc - bmi_mid) / s_bmi))
}

# Transition probabilities for control group
p_1 <- 0.0546     # Healthy to hypertension (control)
p_2 <- rate_to_prob(508.4, 100000, 1)   # Hypertension to CVD (508.4 per 100,000 person-years)
p_3 <- rate_to_prob(20.9, 100, 1)     # CVD to death (20.9 per 100 person-years)

# Intervention case scenario
hr_hta <- 0.45    # Hazard ratio for hypertension with intervention

# Costs
cost_healthy <- 0
cost_hta <- 588.64 # Cost of hypertension management (2024)
cost_stroke <- 6060.10 # Average cost of stroke event (2024)
cost_ami <- 5876.15 # Average cost of AMI event (2024)
cost_cvd <- (cost_stroke + cost_ami) / 2 # Average cost of CVD events
cost_salt <- 3.54 # Cost of salt reduction intervention (2024)
cost_social <- 2.63 # Cost of social intervention (2024)
cost_intervention <- cost_salt + cost_social
cost_death <- 0

# Cycle correction
cycle_length <- 1  # 1 year

# Uptake
uptake <- 0.335  # 33.5% uptake of the intervention

# Disability weights
dw_stroke_mild <- 0.019      # Mean dw for mild stroke
dw_stroke_moderate <- 0.07  # Mean dw for moderate stroke
dw_stroke_severe <- 0.552  # Mean dw for severe stroke
dw_ami <- 0.432 # Mean dw for acute myocardial infarction (AMI) during 1-2 hours
dw_average <- 0.2682233 # See parameters for calculation

# Discount rate
discount_rate <- 0.03

# Load Peru WHO life expectancy table
life_expectancy_table <- load_peru_who_life_table()

# Function to fetch remaining life expectancy for a given age
get_remaining_life <- function(age) {
  age_int <- floor(age)
  # assume table has columns 'age' and 'ex'
  idx <- which(life_expectancy_table$age == age_int)
  if (length(idx) == 0) {
    idx <- which.min(abs(life_expectancy_table$age - age_int))
  }
  life_expectancy_table$ex[idx]
}

# -----------------------------------------------------------
# Create transition matrices
# -----------------------------------------------------------
create_transition_matrix <- function(p_1, hr_hta = 1, p_2_local = p_2, p_3_local = p_3, p_bg) {
  # Four‑state model: Healthy, Hypertension, CVD, Death
  p_1_adj <- 1 - exp(-(p_1 * hr_hta))   # Adjusted probability Healthy → HTA
  matrix(c(
    # From Healthy
    1 - (p_1_adj + p_bg), p_1_adj, 0, p_bg,
    # From Hypertension
    0, 1 - (p_2_local + p_bg), p_2_local, p_bg,
    # From CVD
    0, 0, 1 - p_3_local, p_3_local,
    # From Death (absorbing state)
    0, 0, 0, 1
  ), nrow = 4, byrow = TRUE)
}
# Helper function to get annual background mortality probability from life table,
# subtracting estimated annual CVD mortality rate.
get_background_mortality <- function(age) {
  # This assumes that the life table provides remaining life expectancy (ex) for the given age
  age_int <- floor(age)
  idx <- which(life_expectancy_table$age == age_int)
  if (length(idx) == 0) {
    idx <- which.min(abs(life_expectancy_table$age - age_int))
  }
  ex <- life_expectancy_table$ex[idx]
  mort_allcause <- 1 / ex
  mort_cvd <- p_3 # Approximate annual CVD mortality rate (replace with real data as available)
  p_bg <- max(0, mort_allcause - mort_cvd)
  return(p_bg)
}

# Create matrices
# Logistic penalization of hypertension HR by SBP and DBP
RRR_base   <- 1 - hr_hta
L_bp <- 0.1
sbp_mid <- 140
s_sbp <- 5
dbp_mid <- 90
s_dbp <- 5
sbp_ref_floor <- 100
dbp_ref_floor <- 60
adjust_sbp <- function(sbp) pmax(sbp, sbp_ref_floor)
adjust_dbp <- function(dbp) pmax(dbp, dbp_ref_floor)
adj_sbp <- adjust_sbp(sbp_mean)
adj_dbp <- adjust_dbp(dbp_mean)
frac_sbp <- L_bp + (1 - L_bp) / (1 + exp((adj_sbp - sbp_mid) / s_sbp))
frac_dbp <- L_bp + (1 - L_bp) / (1 + exp((adj_dbp - dbp_mid) / s_dbp))
# Preserve BMI penalization alongside BP
frac_bmi <- frac_bmi_retained(bmi)
f_total  <- frac_sbp * frac_dbp * frac_bmi
hr_hta_eff <- 1 - RRR_base * f_total
hr_mix     <- (1 - uptake) + uptake * hr_hta_eff
# Compute background mortality probability at initial age
p_bg <- 1 / get_remaining_life(age_init)
mat_control      <- create_transition_matrix(p_1, hr_hta = 1, p_bg = p_bg)
mat_intervention <- create_transition_matrix(p_1, hr_hta = hr_mix, p_bg = p_bg)

# -----------------------------------------------------------
# Discounting function
# -----------------------------------------------------------
apply_discount <- function(value, rate, time) {
  value / (1 + rate)^time
}

# -----------------------------------------------------------
# Run deterministic model
# -----------------------------------------------------------
run_deterministic_model <- function(transition_matrix, costs, disability_weights,
                                    cycles, init_state, is_intervention = FALSE,
                                    discount_rate = 0.03) {

  if (is_intervention) {
    base_costs <- c(costs$intervention_healthy,
                    costs$intervention_hta,
                    costs$intervention_cvd,
                    costs$intervention_death)
  } else {
    base_costs <- c(costs$healthy,
                    costs$hta,
                    costs$cvd,
                    costs$death)
  }

  dw_vector <- c(0,
                 disability_weights$hta,
                 disability_weights$cvd,
                 0)

  state_dist <- matrix(0, nrow = cycles + 1, ncol = 4)
  state_dist[1, ] <- init_state

  costs_by_cycle <- numeric(cycles)
  yld_by_cycle   <- numeric(cycles)
  yll_by_cycle   <- numeric(cycles)
  new_hta_by_cycle <- numeric(cycles)
  new_cvd_by_cycle <- numeric(cycles)

  for (t in 1:cycles) {
    if (t > 1) state_dist[t, ] <- state_dist[t - 1, ] %*% transition_matrix

    # --- YLL calculation ---
    prev_state <- if (t == 1) init_state else state_dist[t - 1, ]
    # deaths this cycle = sum over states of transitions to death (column 4)
    deaths <- sum(prev_state * transition_matrix[, 4])
    # Incident cases this cycle
    new_hta_by_cycle[t] <- prev_state[1] * transition_matrix[1, 2]
    new_cvd_by_cycle[t] <- prev_state[2] * transition_matrix[2, 3]
    age_at_cycle <- floor(age_init + t - 1)
    remaining_life <- get_remaining_life(age_at_cycle)
    # discounted YLL for this cycle
    yll_by_cycle[t] <- apply_discount(deaths * remaining_life, discount_rate, t - 1)
    # -----------------------

    cycle_costs_raw <- sum(state_dist[t, ] * base_costs)
    cycle_yld_raw   <- sum(state_dist[t, ] * dw_vector)
    cycle_costs     <- apply_discount(cycle_costs_raw, discount_rate, t - 1)
    yld_by_cycle[t] <- apply_discount(cycle_yld_raw, discount_rate, t - 1)

    costs_by_cycle[t] <- cycle_costs
  }

  total_hta_cases <- sum(new_hta_by_cycle)
  total_cvd_cases <- sum(new_cvd_by_cycle)

  total_costs <- sum(costs_by_cycle)
  total_yld   <- sum(yld_by_cycle)
  total_yll   <- sum(yll_by_cycle)
  total_dalys <- total_yld + total_yll

  list(
    state_dist     = state_dist,
    costs_by_cycle = costs_by_cycle,
    yld_by_cycle   = yld_by_cycle,
    yll_by_cycle   = yll_by_cycle,
    total_costs    = total_costs,
    total_yll      = total_yll,
    total_dalys    = total_dalys,
    new_hta_by_cycle   = new_hta_by_cycle,
    new_cvd_by_cycle   = new_cvd_by_cycle,
    total_hta_cases    = total_hta_cases,
    total_cvd_cases    = total_cvd_cases
  )
}

# -----------------------------------------------------------
# Run both strategies
# -----------------------------------------------------------
costs <- list(
  healthy               = cost_healthy,
  hta                   = cost_hta,
  cvd                   = cost_cvd,
  death                 = cost_death,
  intervention_healthy  = cost_healthy + cost_intervention,
  intervention_hta      = cost_hta + cost_intervention,
  intervention_cvd      = cost_cvd + cost_intervention,
  intervention_death    = cost_death + cost_intervention
)

disability_weights <- list(
  hta = 0,  # Average disability weight for hypertension
  cvd = dw_average
)

init_states <- c(1, 0, 0, 0)  # All start in Healthy state
cycles <- 40  # Lifetime horizon

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

# Calculate averted cases
hta_averted <- control_results$total_hta_cases - intervention_results$total_hta_cases
cvd_averted <- control_results$total_cvd_cases - intervention_results$total_cvd_cases

# Express averted cases as rates

hta_averted_per100   <- hta_averted * 100
cvd_averted_per100  <- cvd_averted * 100

# Compute percent of cases averted
hta_percent_averted <- 100 * hta_averted / control_results$total_hta_cases
cvd_percent_averted <- 100 * cvd_averted / control_results$total_cvd_cases

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
  ICER               = c(NA, icer),
  Averted_HTA        = c(NA, hta_averted),
  Averted_CVD        = c(NA, cvd_averted),
  Averted_HTA_per100   = c(NA, hta_averted_per100),
  Averted_CVD_per100  = c(NA, cvd_averted_per100),
  Percent_HTA_Averted  = c(NA, hta_percent_averted),
  Percent_CVD_Averted  = c(NA, cvd_percent_averted)
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

  # Sample HTA cost and disability weight
  cost_hta_i <- rgamma(1, shape = 476, rate = 1)
  dw_average_i   <- rbeta(1, shape1 = 26.45663, shape2 = 72.17015)

  # Sample HTA→CVD and CVD→Death via gamma distributions (dummy)
  rate_p2_i <- rgamma(1, shape = 14.8, rate = 0.029)
  p2_i <- rate_to_prob(rate_p2_i, 100000, 1) # nolint
  rate_p3_i <- rgamma(1, shape = 37.1, rate = 1.78)
  p3_i <- rate_to_prob(rate_p3_i, 100, 1) # nolint

  # Build cost list (four states)
  costs_i <- list(
    healthy               = cost_healthy,
    hta                   = cost_hta_i,
    cvd                   = cost_cvd,
    death                 = cost_death,
    intervention_healthy  = cost_healthy + cost_intervention,
    intervention_hta      = cost_hta_i + cost_intervention,
    intervention_cvd      = cost_cvd + cost_intervention,
    intervention_death    = cost_death + cost_intervention
  )
  # Build disability weights list
  disability_weights_i <- list(
    hta = 0,
    cvd = dw_average_i
  )

  # Logistic penalization of hypertension HR by SBP and DBP in PSA
  RRR_base_i   <- 1 - hr_hta_i
  L_bp <- 0.1
  sbp_mid <- 140
  s_sbp <- 5
  dbp_mid <- 90
  s_dbp <- 5
  sbp_ref_floor <- 100
  dbp_ref_floor <- 60
  adjust_sbp <- function(sbp) pmax(sbp, sbp_ref_floor)
  adjust_dbp <- function(dbp) pmax(dbp, dbp_ref_floor)
  adj_sbp_i <- adjust_sbp(sbp_mean_i)
  adj_dbp_i <- adjust_dbp(dbp_mean_i)
  frac_sbp_i <- L_bp + (1 - L_bp) / (1 + exp((adj_sbp_i - sbp_mid) / s_sbp))
  frac_dbp_i <- L_bp + (1 - L_bp) / (1 + exp((adj_dbp_i - dbp_mid) / s_dbp))
  # Preserve BMI penalization alongside BP
  frac_bmi_i <- frac_bmi_retained(bmi_i)
  f_total_i  <- frac_sbp_i * frac_dbp_i * frac_bmi_i
  hr_hta_eff_i <- 1 - RRR_base_i * f_total_i
  hr_mix_i     <- (1 - uptake) + uptake * hr_hta_eff_i

  # Get background mortality for this sampled initial age
  p_bg_i <- get_background_mortality(age_init_i)
  # Create transition matrices
  mat_control_i      <- create_transition_matrix(p_1, hr_hta = 1, p_2_local = p2_i, p_3_local = p3_i, p_bg = p_bg_i)
  mat_intervention_i <- create_transition_matrix(p_1, hr_hta = hr_mix_i, p_2_local = p2_i, p_3_local = p3_i, p_bg = p_bg_i)

  # Run deterministic model
  res_control_i <- run_deterministic_model(mat_control_i,
                                           costs_i,
                                           disability_weights_i,
                                           cycles,
                                           init_states,
                                           is_intervention = FALSE,
                                           discount_rate)
  res_int_i     <- run_deterministic_model(mat_intervention_i,
                                           costs_i,
                                           disability_weights_i,
                                           cycles,
                                           init_states,
                                           is_intervention = TRUE,
                                           discount_rate)

  # Incremental cost and DALYs averted
  inc_cost <- res_int_i$total_costs - res_control_i$total_costs
  inc_eff  <- res_control_i$total_dalys - res_int_i$total_dalys

  c(inc_cost = inc_cost, inc_eff = inc_eff)
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

# Peruvian DALY threshold based on GDP per capita (I$) threeshold based in opportunity cost
wtp_threshold_i <- 2128 # Inferior threshold
wtp_threshold_s <- 8372.5 # Superior threshold
# Cost-effectiveness plane
ce_plane_df <- psa_df
# Calculate centroid
centroid_df <- data.frame(
  inc_eff = mean(ce_plane_df$inc_eff, na.rm = TRUE),
  inc_cost = mean(ce_plane_df$inc_cost, na.rm = TRUE)
)
gg_ce_plane <- ggplot(ce_plane_df, aes(x = inc_eff, y = inc_cost)) +
  geom_point(color = "#b20000a3", alpha = 0.4) +
  geom_abline(intercept = 0, slope = wtp_threshold_i, color = "#0d009e", linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 0.6) +
  geom_hline(yintercept = 0, color = "black", size = 0.6) +
  stat_ellipse(type = "norm", level = 0.95, color = "#999999", linetype = "solid", size = 1) +
  geom_point(data = centroid_df, aes(x = inc_eff, y = inc_cost), color = "#999999", size = 4, shape = 8) +
  labs(
    x = "Incremental DALYs Averted",
    y = "Incremental Costs"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
print(gg_ce_plane)

ggsave(gg_ce_plane,
filename = "Cost effectiveness plane.jpeg",
width = 20, height = 20, units = "cm",
dpi = 1600,
path = "output/figs")

# Cost-effectiveness acceptability curve (CEAC)
wtp_vals <- seq(0, 10000, by = 100)
ceac <- sapply(wtp_vals, function(wtp) {
  mean(psa_df$inc_eff - psa_df$inc_cost / wtp >= 0, na.rm = TRUE)
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
  geom_line(color = "#0072B2", size = 1) +
  labs(
    x     = "Willingness-to-pay Threshold",
    y     = "EVPI",
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(gg_evpi)