# ---
# Supplementary materials
# Author: Darwin Del Castillo
# ---

## Creating a Markov Transition Model
# Loading packages
pacman::p_load(DiagrammeR,
               DiagrammeRsvg,
               rsvg)

# Creating the main analysis graph
library(DiagrammeR)

myGraphSimple <- grViz("
digraph MarkovModel {
  rankdir=LR;

  ## Set default state style
  node [shape = circle, style=solid, 
  fontname='Helvetica', fontsize=8, 
  width=0.6, height=0.6];

  ## Define states
  Healthy           [label = 'Healthy'];
  Hypertension      [label = 'Hypertension'];
  CVD               [label = 'CV events'];
  Death             [label = 'Death']
  
  ## Edge defaults
  edge [fontname='Helvetica', fontsize=8];

  ## Define transitions with labels
  Healthy -> Healthy           [label = 'C'];
  Healthy -> Hypertension      [label = 'p_1'];
  Healthy -> Death             [label = 'p_2'];

  Hypertension -> Hypertension [label = 'C'];
  Hypertension -> CVD          [label = 'p_3'];
  Hypertension -> Death        [label = 'p_4'];
  
  CVD -> Death [label = 'p_5'];
  CVD -> CVD [label = 'C'];

  
  Death -> Death [label = '1'];
}
")

myGraph <- grViz("
digraph MarkovModel {
  rankdir=LR;

  ## Set default state style
  node [shape = circle, style=solid, 
  fontname='Helvetica', fontsize=8, 
  width=0.6, height=0.6];

  ## Define states
  Healthy           [label = 'Healthy'];
  Prehypertension   [label = 'Prehypertension'];
  Hypertension      [label = 'Hypertension'];
  CVD               [label = 'CV events'];
  Death             [label = 'Death']
  
  ## Edge defaults
  edge [fontname='Helvetica', fontsize=8];

  ## Define transitions with labels
  Healthy -> Healthy           [label = 'C'];
  Healthy -> Prehypertension   [label = 'p_2'];
  Healthy -> Hypertension      [label = 'p_3'];
  Healthy -> CVD               [label = 'p_4'];
  Healthy -> Death             [label = 'p_5'];

  Prehypertension -> Healthy                  [label = 'p_6'];
  Prehypertension -> Hypertension             [label = 'p_7'];
  Prehypertension -> Prehypertension          [label = 'C'];
  Prehypertension -> CVD                      [label = 'p_9'];

  Hypertension -> Hypertension [label = 'C'];
  Hypertension -> CVD          [label = 'p_10'];
  Hypertension -> Death        [label = 'p_11'];
  
  CVD -> Death [label = 'p_15'];
  
  Death -> Death [label = '1'];
}
")

myGraphSex <- grViz("
digraph MarkovSexSplit {
  rankdir = LR;

  # Shared styling
  node [shape = circle, fontname = Helvetica, 
  fontsize = 8, width = 0.6, height = 0.6];
  edge [fontname = Helvetica, fontsize = 8];

  ############################
  ## Male sub-graph
  subgraph cluster_m {
    label = 'Male (M)';
    style = dashed;

    Healthy_M         [label = 'Healthy'];
    PreHT_M           [label = 'Prehypertension'];
    HT_M              [label = 'Hypertension'];
    CVD_M             [label = 'CV events'];
    HistCVD_M         [label = 'History of CV events'];
    Death_M           [label = 'Death'];

    Healthy_M -> Healthy_M       [label = 'p1_M'];
    Healthy_M -> PreHT_M         [label = 'p2_M'];
    Healthy_M -> HT_M            [label = 'p3_M'];
    Healthy_M -> CVD_M           [label = 'p4_M'];
    Healthy_M -> Death_M         [label = 'p5_M'];

    PreHT_M -> HT_M              [label = 'p6_M'];
    PreHT_M -> PreHT_M           [label = 'p7_M'];
    PreHT_M -> CVD_M             [label = 'p8_M'];

    HT_M -> HT_M                 [label = 'p9_M'];
    HT_M  -> CVD_M               [label = 'p10_M'];
    HT_M  -> Death_M             [label = 'p11_M'];

    CVD_M  -> HistCVD_M          [label = 'p12_M'];
    HistCVD_M -> HistCVD_M       [label = 'p13_M'];
    HistCVD_M -> Death_M         [label = 'p14_M'];
    CVD_M  -> Death_M            [label = 'p15_M'];

    Death_M -> Death_M           [label = '1'];
  }

  ############################
  ## Female sub-graph
  subgraph cluster_f {
    label = 'Female (F)';
    style = dashed;

    Healthy_F         [label = 'Healthy'];
    PreHT_F           [label = 'Prehypertension'];
    HT_F              [label = 'Hypertension'];
    CVD_F             [label = 'CV events']
    Death_F           [label = 'Death'];

    Healthy_F -> Healthy_F       [label = 'p1_F'];
    Healthy_F -> PreHT_F         [label = 'p2_F'];
    Healthy_F -> HT_F            [label = 'p3_F'];
    Healthy_F -> CVD_F           [label = 'p4_F'];
    Healthy_F -> Death_F         [label = 'p5_F'];

    PreHT_F -> HT_F              [label = 'p6_F'];
    PreHT_F -> PreHT_F           [label = 'p7_F'];
    PreHT_F -> CVD_F             [label = 'p8_F'];

    HT_F -> HT_F                 [label = 'p9_F'];
    HT_F  -> CVD_F               [label = 'p10_F'];
    HT_F  -> Death_F             [label = 'p11_F'];

    CVD_F  -> Death_F            [label = 'p15_F'];

    Death_F -> Death_F           [label = '1'];
  }
}
")

print(myGraphSimple)

# Converting the graph into a JPEG file
graph_svg <- export_svg(myGraphSimple)

rsvg_png(charToRaw(graph_svg),
         file = "output/figs/Simple Markov Model.png",
         width = 1400,
         height = 800)

# Creating the sensitivity analyses graphs

# Creating the first sensitivity analysis graph including inpatient costs
myGraph_S1 <- grViz("
digraph MarkovModelS1 {
  rankdir=LR;

  ## Set default state style
  node [shape = circle, style=solid, 
  fontname='Helvetica', fontsize=8, 
  width=0.6, height=0.6];

  ## Define states
  Healthy           [label = 'Healthy'];
  Prehypertension   [label = 'Prehypertension'];
  Hypertension      [label = 'Hypertension'];
  CVD               [label = 'CV events'];
  Death             [label = 'Death'];
  Hospitalization   [label = 'Hospitalization']
  
  ## Edge defaults
  edge [fontname='Helvetica', fontsize=8];

  ## Define transitions with labels
  Healthy -> Healthy           [label = 'p_1'];
  Healthy -> Prehypertension   [label = 'p_2'];
  Healthy -> Hypertension      [label = 'p_5'];
  Healthy -> CVD               [label = 'p_6'];
  Healthy -> Death             [label = 'p_7'];

  Prehypertension -> Hypertension             [label = 'p_4'];
  Prehypertension -> Prehypertension          [label = 'p_3'];
  Prehypertension -> CVD                      [label = 'p_8'];
  Prehypertension -> Death                    [label = 'p_9'];
  
  Hypertension -> Hypertension [label = 'p_10'];
  Hypertension -> CVD          [label = 'p_11'];
  Hypertension -> Death        [label = 'p_12'];
  Hypertension -> Hospitalization [label = 'p_12_1'];
  
  CVD -> CVD [label = 'p_13'];
  CVD -> Hospitalization [label = 'p_13_1'];
  CVD -> Death [label = 'p_14'];
  
  Hospitalization -> CVD [label = 'p_13_2'];
  Hospitalization -> Death [label = 'p_16'];
  Hospitalization -> Hypertension [label = 'p_12_2'];
  
  Death -> Death [label = '1'];
}
")

print(myGraph_S1)

# Converting the graph into a JPEG file
graph_svgS1 <- export_svg(myGraph_S1)

rsvg_png(charToRaw(graph_svgS1),
         file = "output/figs/Markov Model Graph S1.png",
         width = 2000,
         height = 1000)

# Creating the second sensitivity analysis graph including outpatient costs
myGraph_S2 <- grViz("
digraph MarkovModelS2 {
  rankdir=LR;

  ## Set default state style
  node [shape = circle, style=solid, 
  fontname='Helvetica', fontsize=8, 
  width=0.6, height=0.6];

  ## Define states
  Healthy           [label = 'Healthy'];
  Prehypertension   [label = 'Prehypertension'];
  Hypertension      [label = 'Hypertension'];
  CVD               [label = 'CV events'];
  Death             [label = 'Death'];
  Outpatient_Costs         [label = 'Medicines']
  
  ## Edge defaults
  edge [fontname='Helvetica', fontsize=8];

  ## Define transitions with labels
  Healthy -> Healthy           [label = 'p_1'];
  Healthy -> Prehypertension   [label = 'p_2'];
  Healthy -> Hypertension      [label = 'p_5'];
  Healthy -> CVD               [label = 'p_6'];
  Healthy -> Death             [label = 'p_7'];

  Prehypertension -> Hypertension             [label = 'p_4'];
  Prehypertension -> Prehypertension          [label = 'p_3'];
  Prehypertension -> CVD                      [label = 'p_8'];
  Prehypertension -> Death                    [label = 'p_9'];
  
  Hypertension -> Hypertension [label = 'p_10'];
  Hypertension -> CVD          [label = 'p_11'];
  Hypertension -> Death        [label = 'p_12'];
  Hypertension -> Outpatient_Costs [label = 'p_12_1'];
  
  CVD -> CVD [label = 'p_13'];
  CVD -> Outpatient_Costs [label = 'p_13_1'];
  CVD -> Death [label = 'p_14'];
  
  Outpatient_Costs -> Outpatient_Costs [label = 'p_13_2'];
  Outpatient_Costs -> Death [label = 'p_16'];

  Death -> Death [label = '1'];
}
")

print(myGraph_S2)

graph_svgS2 <- export_svg(myGraph_S2)


rsvg_png(charToRaw(graph_svgS2),
         file = "output/figs/Markov Model Graph S2.png",
         width = 2000,
         height = 1000)