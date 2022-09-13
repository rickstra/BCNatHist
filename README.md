
# BCNatHist: A class of continuous tumor growth natural history models of breast cancer. 

## Description

[[1]](#1)

### To be added:

  * Implemented functions for various model risk predictions.
  * Microsimulation framework. 
  * Metastasis extension. 

## Installation

devtools::install_github()

### Dependencies

foreach, doParallel, itertools.

## How-to

Start by defining a natural history model by calling `DefineBCModel`. Each of the four submodels (`onset, growth, sympt, sens`) can be specified with a formula, e.g. `onset = ~ x1 + x2`. Make sure that the variable names exist in the intended data set.  

  > `model <- DefineBCModel(onset = ~ x1 + x2, sens = ~ x3)`

Once defined, the model can be fitted to the `data` using maximum likelihood by calling `EstimateBCModel`. The function returns the model object updated with the fitted parameter values and the estimated covariance matrix (by default).  

  > `model <- EstimateBCModel(model, data)`

The estimated parameter values can then be summarized, including Wald confidence intervals based on the estimated covariance matrix.  

  > `summary(model)`

If the parameter values are already known or otherwise to be fixed as `par_values`, they can be manually added to the model by overwriting the `model$par` variable, or added when defining the model using `DefineBCmodel`. When overriding parameter values, also change the `model$fitted` variable to `TRUE`.   

  > `model$par <- par_values`  
  > `model$fitted <- TRUE`  
  > \# or  
  > `model <- DefineBCModel(start.par = par_values)`  
  > `model$fitted <- TRUE`  

## References

<a id="1">[1]</a> 
J. R. Strandberg, K. Humphreys.  
Statistical models of tumour onset and growth for modern breast cancer screening cohorts.  
Mathematical Biosciences 318 (2019) 108270.  
[doi:10.1016/j.mbs.2019.108270](https://doi.org/10.1016/j.mbs.2019.108270).

<a id="2">[2]</a> 
Strandberg R, Czene K, Eriksson M, Hall P, Humphreys K.  
Estimating Distributions of Breast Cancer Onset and Growth in a Swedish Mammography Screening Cohort.  
Cancer Epidemiol Biomarkers Prev (2022) 31(3):569–577.  
[doi: 10.1158/1055-9965.EPI-21-1011](https://doi.org/10.1158/1055-9965.EPI-21-1011).
