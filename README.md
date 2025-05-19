# Cost-effectiveness analysis of a salt substitute intervention in Northern Peru
This repository analyzes the economic impact of a salt substitute intervention in northern Peru. The data belongs to the CRONICAS Center of Excellence in Chronic Diseases at the Universidad Peruana Cayetano Heredia (and the study participants) and is available upon reasonable request to its owners.

## Packages used
- tidyverse
- haven
- heemod
- DiagrammeR
- skimr
- labelled
- data.table
- openxlsx
- lme4
- survival
- rigr

## To-do list
### Primary data
- [ ] Compare dichotomic outcomes with sbp and dbp measurements in each round.
- [ ] Estimate the effect of the intervention oon the incidence of pre-HTA.
- [ ] Estimate the effect of the intervention on the incidence of HTA among people with pre-HTA at baseline.
- [ ] Retrieve better estimates of mortality (prob. using SINADEF access)
- [ ] Improve estimates for CVD events. 
- [ ] Introduce the impact of BMI, SBP, and DBP on the HR estimate.

### Markov model
- [ ] Run the model with more flexible setting-up (e.g. rJAGS) 
- [ ] Run PSA analyses, propagating uncertainty with 10000 iterations
- [ ] Improve the model capturing other potential states (e.g. pre-HTA)
- [ ] Capture the independent effect of each CV event (e.g. stroke, MI, etc.)
- [ ] Capture the effect of the intervention across women and men
- [ ] Improve the code to consider age initial as discrete.