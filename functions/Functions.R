# ---
# Functions
# AUTHOR: Darwin Del Castillo
# Date: `r Sys.Date()`
# ---

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
  
  # Convert numeric vector into data frame with age and life expectancy
  life_expectancy_df <- data.frame(
    age = 0:max_age,
    ex  = life_expectancy,
    stringsAsFactors = FALSE
  )
  return(life_expectancy_df)
}

# ------------------------------------------------------------
# Function to convert rate to probability
# ------------------------------------------------------------

rate_to_prob <- function(r, per, to) {
  1 - exp(-(r / per) * to)
}

# ------------------------------------------------------------
# Penalization functions
# ------------------------------------------------------------

# TBD