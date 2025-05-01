```{r additional settings}
# setting the number of iterations
sim_run <- 10000
# disable scientific notation
options(scipen = 999999,
        max.print = 99999999)
```

```{r loading packages}
pacman::p_load(tidyverse,
               haven,
               heemod,
               skimr,
               labelled,
               data.table,
               openxlsx,
               lme4,
               nlme,
               survival,
               rigr)

# external functions
source("functions/Functions.R")
```

```{r importing datasets}
# raw long data
raw_long_dta <- read_dta("data-raw/long.dta") |> 
  as.data.table() |> 
  mutate(across(where(is.labelled), ~ as_factor(.x, 
                                                levels = "labels")),
         across(everything(), ~ {attr(.x, 
                                      "format.stata") <- NULL
         .x})) |>
  mutate(across(everything(), ~`attr<-`(.x, "label", NULL)))

# raw wide data
raw_wide_dta <- read_dta("data-raw/wide.dta") |> 
  as.data.table() |>
  mutate(across(where(is.labelled), ~ as_factor(.x, 
                                                levels = "labels")),
         across(everything(), ~ {attr(.x, 
                                      "format.stata") <- NULL
         .x})) |>
  mutate(across(everything(), ~`attr<-`(.x, "label", NULL)))
```

```{r exploring databases}
# database structure
str(raw_long_dta)
str(raw_wide_dta)

# database names imported correctly
names(raw_long_dta)
names(raw_wide_dta)

# skimming databases
skim(raw_long_dta)
skim(raw_wide_dta)
```

```{r describing data}
# desenlaces: variables con x son autorreporte y las que no tienen x son medidas en el estudio
## HTA
table(raw_long_dta$ht)
table(raw_long_dta$htdx) 
table(raw_long_dta$ht5) 
table(raw_long_dta$htdxtx) 

## Pre-HTA
table(raw_long_dta$preht)
table(raw_long_dta$preht5)
```

```{r replicating analysis, include = FALSE, evaluate = FALSE}
## survival
### limpieza de fechas
# raw_wide_dta <- raw_wide_dta |>
#  mutate(fecha_inicio = as.Date("2014-08-08"), #fecha de inicio de intervencion: 8 agosto 2014
#         fecha_final = as.Date(f6fechaent, origin = "1970-01-01"), # convertir f6fechaent a Date
#         dtime = as.numeric(fecha_final - fecha_inicio), # calcular tiempo en días
#         dtime = ifelse(dtime < 0, 0, dtime))

# nohypertension_baseline <- subset(raw_wide_dta, ht5 == "No")
# Removing NA from f6ht5
# nohypertension_baseline <- nohypertension_baseline[!is.na(nohypertension_baseline$f6ht5) &                                             !is.na(nohypertension_baseline$dtime), ]
# cox_model <- coxph(Surv(dtime, f6ht5) ~ factor(intervencion) + 
#                     frailty(codigovilla), 
#                     data = nohypertension_baseline,
#                     na.action = na.exclude)
#### cleaning wide data for survival analysis
```
