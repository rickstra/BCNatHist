#' @title DefineBCModel
#' @description Set up a natural history model.
#' @param onset  Specifies covariates to include in the onset submodel.
#' @param growth Specifies covariates to include in the inverse growth rate submodel.
#' @param sympt  Specifies covariates to include in the symptomatic detection submodel.
#' @param sens   Specifies covariates to include in the screening sensitivity submodel.
#' @param onset_dep_growth Indicate if the inverse growth rate should depend on the age at onset (Strandberg et al. 2022).
#' @param data An optional dataframe representing the data the model is to be fitted to. Necessary if e.g. the formulas include spline functions. If provided, all predictor names in the four submodels must exist. 
#' @param start_par Optionally specify (starting) parameter values.
#' @return A list-type object of class \code{"BCmodel"}.
#' \describe{
#'   \item{$onset}{The specified formula for the onset submodel.}
#'   \item{$growth}{The specified formula for the growth rate submodel.}
#'   \item{$sympt}{The specified formula for the symptomatic detection submodel.}
#'   \item{$sens}{The specified formula for the screening sensitivity submodel.}
#'   \item{$onset_dep_growth}{Indicates if the mean of the inverse growth rate distribution depends on the age at onset}
#'   \item{$par}{Numeric vector containing parameter values for all four submodels, either provided by start.par or assigned by the function.}
#'   \item{$fitted}{Indicates if the model parameters have been fitted to data. Set to FALSE.}
#' }
#' @examples
#' DefineBCModel()
#' DefineBCModel(sens = ~ percent_density, onset_dep_growth = TRUE)
#' DefineBCModel(onset = ~ family_history + birads, growth = ~ bmi)
#' @export
DefineBCModel <- function(onset = ~ 1, growth = ~ 1, sympt = ~ 1, sens = ~ 1, 
                          onset_dep_growth = FALSE, data = NULL, start_par = NULL) {
  
  if (is.null(data)) {
    var_names <- c(all.vars(onset), all.vars(growth), all.vars(sympt), all.vars(sens))
    data <- as.data.frame(matrix(1, 1, length(var_names)))
    names(data) <- var_names
  }
  
  o_names <- colnames(model.matrix(onset, data))[-1]
  g_names <- colnames(model.matrix(growth, data))[-1]
  sym_names <- colnames(model.matrix(sympt, data))[-1]
  sen_names <- colnames(model.matrix(sens, data))[-1]
  
  par_names <- c("A", "B", "phi", "beta1", "delta0", o_names, "mu0", g_names, "eta0", sym_names, "beta0", sen_names)
  
  if(is.null(start_par)){
    start_par <- c(-2.28, -7.09, -0.49, 0.47,
                   -3.08, sample(c(-0.01, 0.01), length(o_names), T),
                   0.07, sample(c(-0.01, 0.01), length(g_names), T),
                   -8.85, sample(c(-0.01, 0.01), length(sym_names), T),
                   -5.37, sample(c(-0.01, 0.01), length(sen_names), T))
  }
  if (onset_dep_growth) {
    par_names <- c(par_names, "odg")
    start_par <- c(start_par, 0.01)
  } 
  names(start_par) <- par_names
  model <- list(onset = onset, growth = growth, sympt = sympt, sens = sens, 
                onset_dep_growth = onset_dep_growth, 
                par = start_par, fitted = FALSE)
  class(model) <- "BCModel"
  return(model)
}
