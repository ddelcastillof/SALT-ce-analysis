# ---
# Additional
# Author: Darwin Del Castillo
# ---

## Creating a Markov Transition Model
# Loading packages
pacman::p_load("DiagrammeR",
               "DiagrammeRsvg",
               "rsvg")

# Creating the graph
myGraph <- grViz("
digraph MarkovModel {
  rankdir=LR;

  // Set default state style
  node [shape = circle, style=solid, fontname='Helvetica', fontsize=8, width=0.6, height=0.6];

  // Define states
  Healthy       [label = 'Healthy'];
  Hypertension  [label = 'Hypertension'];
  CVEvent       [label = 'CV Event'];
  Dead          [label = 'Dead'];
  
  // Edge defaults
  edge [fontname='Helvetica', fontsize=8];

  // Define transitions with labels
  Healthy -> Healthy        [label = 'p_1'];
  Healthy -> Hypertension   [label = 'p_2'];
  Healthy -> Dead           [label = 'p_3'];

  Hypertension -> Hypertension [label = 'p_4'];
  Hypertension -> Dead         [label = 'p_5'];
  Hypertension -> CVEvent      [label = 'p_6'];
  
  CVEvent -> Dead              [label = 'p_7'];
  CVEvent -> CVEvent           [label = 'p_8'];

  Dead -> Dead                 [label = '1'];
}
")

print(myGraph)

## Converting the graph into a SVG file
graph_svg <- export_svg(myGraph)

rsvg_png(charToRaw(graph_svg),
         file = "figs/Markov Model Graph.png",
         width = 2000,
         height = 1600)

