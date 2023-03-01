#' @title BuildParameters
#' @description Predict BC risk with an upper limit on tumor size, both if and if not attending the next screening.
#' 
#' @param par Vector of parameter values, transformed to have support on R. 
#' @param model Natural history model of class \code{BCModel} as defined by \code{DefineBCModel}.
#' @param data Dataframe containing the covariate data of the model. One entry per individual.
#'
#' @return \code{matrix} with the following variables matched from the input data:
#' \describe{
#' \item{[, 1]} A
#' \item{[, 2]} B
#' \item{[, 3]} delta
#' \item{[, 4]} phi
#' \item{[, 5]} mu
#' \item{[, 6]} eta
#' \item{[, 7]} beta_s
#' \item{[, 8]} beta
#' \item{[, 9]} onset-dependent growth coefficient (optional)
#' @examples
#' 
#' @export
BuildParameters <- function(par, model, data) {
  
  o_pred <- model.matrix(model$onset, data)
  g_pred <- model.matrix(model$growth, data)
  sym_pred <- model.matrix(model$sympt, data)
  scr_pred <- model.matrix(model$sens, data)
  
  # Core
  A <- -exp(par[1])
  B <- exp(par[2])
  phi <- exp(par[3])
  beta <- par[4]
  par <- par[-(1:4)]
  # Onset
  o_par <- par[1:ncol(o_pred)]
  par <- par[-(1:ncol(o_pred))]
  # Growth
  g_par <- par[1:ncol(g_pred)]
  par <- par[-(1:ncol(g_pred))]
  # Symptomatic detection
  sym_par <- par[1:ncol(sym_pred)]
  par <- par[-(1:ncol(sym_pred))]
  # Screen detection
  scr_par <- par[1:ncol(scr_pred)]
  par <- par[-(1:ncol(scr_pred))]
  
  # linear predictors
  o_pred <- o_pred %*% o_par
  g_pred <- g_pred %*% g_par
  sym_pred <- sym_pred %*% sym_par
  scr_pred <- scr_pred %*% scr_par
  
  par_names <- c("A", "B", "delta", "phi", "mu", "eta", "beta1", "beta2")
  
  # Onset-dependent growth?
  if(model$onset_dep_growth) {
    par_names <- c(par_names, "odg")
  } else {
    par <- NULL
  }
  par_mat <- cbind(A, B, exp(o_pred), phi, exp(g_pred), 
                   exp(sym_pred), beta, scr_pred, par)
  colnames(par_mat) <- par_names
  return(par_mat)
}
