# ---
# Cleaning code
# Author: Darwin Del Castillo
# ---

# --- Load required libraries -------------
pacman::p_load(tidyverse,
               haven,
               skimr,
               labelled,
               data.table,
               openxlsx
               )

# --- Load external functions -------------
source("functions/Functions.R")

# -----------------------------------------------------------
# Load SALT data 
# -----------------------------------------------------------

## Raw wide data
raw_wide_dta <- read_dta("data-raw/wide.dta") |> 
  as.data.table() |>
  mutate(across(where(is.labelled), ~ as_factor(.x, 
                                                levels = "labels")),
         across(everything(), ~ {attr(.x, 
                                      "format.stata") <- NULL
         .x})) |>
  mutate(across(everything(), ~`attr<-`(.x, "label", NULL)))

# -----------------------------------------------------------
# Cleaning wide data 
# -----------------------------------------------------------

### Selecting variables from raw data
clean_wide <- raw_wide_dta |>
  select(#Time_0
         entvilla,
         codigo,
         codigogen,
         codigovilla,
         codigovivienda,
         codigofam,
         codigopersona,
         fechaent,
         preht5,
         ht5,
         intervencion,
         #Time_1
         f1fecha_b,
         f1preht5,
         f1ht5,
         f1intervencion,
         #Time_2
         f2fecha_b,
         f2preht5,
         f2ht5,
         f2intervencion,
         #Time_3
         f3fecha_b,
         f3preht5,
         f3ht5,
         f3c_interv,
         #Time_4
         f4fecha_b,
         f4c_interv,
         f4preht5,
         f4ht5,
         #Time_5
         f5fecha_b,
         f5c_interv,
         f5preht5,
         f5ht5,
         #Time_6
         f6fecha,
         f6c_interv,
         f6preht5,
         f6ht5) |>
mutate(across(c(f1fecha_b, f2fecha_b, f3fecha_b,
                f4fecha_b, f5fecha_b, f6fecha),
              ~ as.Date(.x, format("%d/%m/%Y")))) |>
mutate(f1intervencion = factor(f1intervencion,
                        levels = c(0, 1),
                        labels = c("Villa NO Intervenida",
                                      "Villa Intervenida")),
       f2intervencion = factor(f2intervencion,
                        levels = c(0, 1, 2),
                        labels = c("Villa NO Intervenida",
                                   "Villa Intervenida",
                                   "Villa Intervenida"))
  ) |>
mutate(intervencion = factor(intervencion,
                        levels = c(0),
                        labels = c("Villa NO Intervenida"))) |>
rename(#Time_0
         f0visit_date = fechaent,
         f0preht = preht5,
         f0ht = ht5,
         f0intervention = intervencion,
         #Time_1
         f1visit_date = f1fecha_b,
         f1preht = f1preht5,
         f1ht = f1ht5,
         f1intervention = f1intervencion,
         #Time_2
         f2visit_date = f2fecha_b,
         f2preht = f2preht5,
         f2ht = f2ht5,
         f2intervention = f2intervencion,
         #Time_3
         f3visit_date = f3fecha_b,
         f3preht = f3preht5,
         f3ht = f3ht5,
         f3intervention = f3c_interv,
         #Time_4
         f4visit_date = f4fecha_b,
         f4preht = f4preht5,
         f4ht = f4ht5,
         f4intervention = f4c_interv,
         #Time_5
         f5visit_date = f5fecha_b,
         f5preht = f5preht5,
         f5ht = f5ht5,
         f5intervention = f5c_interv,
         #Time_6
         f6visit_date = f6fecha,
         f6preht = f6preht5,
         f6ht = f6ht5,
         f6intervention = f6c_interv)
# -----------------------------------------------------------
# Cleaning variables
# -----------------------------------------------------------
## preht is coded as factor No, Yes and Pre-ht5, so I will recode it as Yes, No
## This problem persist across all time points
clean_wide <- clean_wide |>
  mutate(across(c(f0preht, f1preht, f2preht, 
                  f3preht, f4preht, f5preht, f6preht),
                ~ fct_recode(.x, "Yes" = "Pre-ht5")))

#-----------------------------------------------------------
# Pre-processing data
#-----------------------------------------------------------
## Reshaping the data to long format

### Identifying variables not to be pivoted
id_cols <- c("entvilla", "codigo", "codigogen",
             "codigovilla", "codigovivienda",
             "codigofam", "codigopersona")
### Pivoting the data to long format
clean_long <- clean_wide |>
  pivot_longer(
    cols = -all_of(id_cols),
    names_to = c("time", ".value"),
    names_pattern = "(f[0-9])?(.+)"
  )

#-----------------------------------------------------------
# Defining time varying exposure
#-----------------------------------------------------------

## Defining time varying exposure
intervention_schedule <- data.frame(
  codigovilla = c("020", "018", "027", "107", "012", "010"),
  switch_date = as.Date(c("2014-08-08",  # First village
                         "2015-02-12",   # Second village
                         "2015-07-14",   # Third village
                         "2016-01-19",   # Fourth village
                         "2016-05-22",   # Fifth village
                         "2016-10-22"))  # Sixth village
)

# Merging with clean_long to obtain the intervention schedule
clean_long <- clean_long |>
  left_join(intervention_schedule, by = "codigovilla")

# Back transforming intervention variable to numeric
clean_long <- clean_long |>
  mutate(
    intervention_numeric = ifelse(intervention == "Villa Intervenida", 1, 0)
  )
# Obtaining time survival variable

clean_long <- clean_long |>
  arrange(codigo, time) |>
  group_by(codigo) |>
  mutate(
    # Transforming time variable to numeric
    time_num = as.numeric(str_remove(time, "f")),
    # Calculating time difference between visits
    time_diff = as.numeric(difftime(visit_date, lag(visit_date, default = first(visit_date)), units = "days")),
    # Handling missing time_diff
    time_diff = ifelse(time_num == 0, 0, time_diff),
    # Calculating cumulative time
    cum_time = cumsum(time_diff),
    # Calculating time since intervention
    tstart = lag(cum_time, default = 0),
    tstop = cum_time
  ) |>
  ungroup()

#-----------------------------------------------------------
# Removing objects no longer needed
#-----------------------------------------------------------
rm(id_cols, intervention_schedule)

#-----------------------------------------------------------
# Sub-setting participants from long format data with pre-hypertension at time 0
#-----------------------------------------------------------
prehypertension_baseline <- clean_wide |>
  filter(!f0preht == "Yes")

survival_data <- clean_wide |>
  filter(f0ht == "No")