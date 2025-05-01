# ---
# Name: Data dictionary
# Author: Darwin Del Castillo
# ---

## Creating dictionary files
# Loading packages
pacman::p_load(tidyverse,
               haven,
               openxlsx,
               labelled,
               data.table)

# Importing datasets
raw_long_dta <- read_dta("data-raw/long.dta") |> 
  as.data.table() |> 
  mutate(across(where(is.labelled), ~ as_factor(.x, 
                                                levels = "labels")),
         across(everything(), ~ {attr(.x, 
                                      "format.stata") <- NULL
         .x}))

raw_wide_dta <- read_dta("data-raw/wide.dta") |> 
  as.data.table() |>
  mutate(across(where(is.labelled), ~ as_factor(.x, 
                                                levels = "labels")),
         across(everything(), ~ {attr(.x, 
                                      "format.stata") <- NULL
         .x}))
data_dictionary_long <- data.frame(
  variable_name = names(raw_long_dta),
  variable_label = sapply(raw_long_dta, function(x) attr(x, "label")),
  row.names = NULL
)

# Adding an additional line because the phq-9 variable does not have a label attribute
data_dictionary_wide <- data.frame(
  variable_name = names(raw_wide_dta),
  variable_label = sapply(raw_wide_dta, function(x) {
    label <- attr(x, "label")
    if (is.null(label)) "" else label
  }),
  row.names = NULL
)

# Saving dictionaries
salt_dictionary <- createWorkbook()
addWorksheet(salt_dictionary, "data_dictionary_long")
addWorksheet(salt_dictionary, "data_dictionary_wide")
writeData(wb = salt_dictionary, data_dictionary_long, sheet = "data_dictionary_long", startRow = 1, startCol = 1)
writeData(wb = salt_dictionary, data_dictionary_wide, sheet = "data_dictionary_wide", startRow = 1, startCol = 1)
saveWorkbook(salt_dictionary, "data-raw/dictionary/salt_dictionary.xlsx", overwrite = TRUE)

# Removing from memory
rm(salt_dictionary, data_dictionary_wide, data_dictionary_long)