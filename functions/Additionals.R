# ---
# Additional
# Author: Darwin Del Castillo
# ---

## Creating a Markov Transition Model
# Loading packages
pacman::p_load(DiagrammeR,
               DiagrammeRsvg,
               rsvg)

# Creating the graph
myGraph <- grViz("
digraph MarkovModel {
  rankdir=LR;

  ## Set default state style
  node [shape = circle, style=solid, fontname='Helvetica', fontsize=8, width=0.6, height=0.6];

  ## Define states
  Healthy           [label = 'Healthy'];
  Prehypertension  [label = 'Prehypertension'];
  Hypertension      [label = 'Hypertension']
  
  ## Edge defaults
  edge [fontname='Helvetica', fontsize=8];

  ## Define transitions with labels
  Healthy -> Healthy           [label = 'p_1'];
  Healthy -> Prehypertension  [label = 'p_2'];
  Healthy -> Hypertension      [label = 'p_3'];

  Prehypertension -> Hypertension             [label = 'p_4'];
  Prehypertension -> Prehypertension         [label = 'p_5'];
  
  Hypertension -> Hypertension [label = '1']
}
")

print(myGraph)

# Converting the graph into a SVG file
graph_svg <- export_svg(myGraph)

rsvg_png(charToRaw(graph_svg),
         file = "figs/Markov Model Graph.png",
         width = 2000,
         height = 1600)

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
write.xlsx(data_dictionary_long, file = "data-raw/dictionary/data_dictionary_long.xlsx")

write.xlsx(data_dictionary_wide, file = "data-raw/dictionary/data_dictionary_wide.xlsx")

# Removing from memory
rm(data_dictionary_long, data_dictionary_wide)