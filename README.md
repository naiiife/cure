# Causal inference for cure models

## curefunctions.R

Estimate survival function in always-uncured
 
surv_uu(): proposed method based on binary substitutional variable
 
surv_uu1(): based on non-binary substitutional variable (not used in our study)
 
surv_pi(): proposed method based on principal ignorability
 
surv_uu_ci(): survival functions with bootstrap confidence intervals
 
newdata = NULL: estimate Suu(t)
 
newdata = c(x,v,w): estimate Suu(t|x,w)

## simulation.R

Simulation studies and sensitivity analysis

## application.R

Real data application, estimating strata proportions, survival functions and sensitivity analysis


  \n


To cite:

Wang Y, Deng Y and Zhou X-H. Causal inference for time-to-event data with a cured subpopulation. _Biometrics_. 2024.
