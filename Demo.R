## Optionally install dependencies
# install.packages(c("foreach", "doParallel", "itertools", "Rcpp", "RcppAmradillo"))

## Choose either a) or b):
## a) Build and install package from source
# devtools::install_github("rickstra/BCNatHist/BCNatHist", build = TRUE)
## b) Install package from pre-built tarball
# install.packages("https://github.com/rickstra/BCNatHist/raw/main/BCNatHist_1.0.tar.gz", repos = NULL)

#-------------------------------------------------------------------------------

## load example data set
load(url("https://github.com/rickstra/BCNatHist/raw/main/demo_data.RData"))

## Define, fit, and summarise an example model

demo_model <- BCNatHist::DefineBCModel(sens = ~pd)
# demo_model$par <- c(-2.4564, -6.7171, -0.5466, 0.4694, -2.8179, 0.0136, -8.8572, -4.9461, -2.1686)
demo_model <- BCNatHist::EstimateBCModel(demo_model, demo_data, n_core, gauss_kronrod_set = 4, gauss_leguerre_set = 2, 
                                         n_cores = parallel::detectCores() - 1) # Mind the number of cores you want to use
summary(demo_model)
