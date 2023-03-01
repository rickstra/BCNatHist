#' @title EstimateBCModel
#' @description Fit a BCModel to cohort data.
#'
#' @param model Natural history model of class \code{BCModel} as defined by \code{DefineBCModel}.
#' @param data Dataframe containing the cohort data. One entry per individual.
#' @param base_variables List linking variable names of \code{data} to variables used in the prediction: 
#' \describe{
#'   \item{entry}{The age at the start of the prediction.}
#'   \item{exit}{The age at end of prediction follow-up.}
#'   \item{is_case}{BC case status (case = 1).}
#'   \item{mode}{Mode of detection (screen-detected = 1).}
#'   \item{size}{tumor size (mm).}
#'   \item{scr_hist}{Screening histories (each entry a list of ages screened).}
#' }
#' @param hessian Indicates if the Hessian matrix of the maximum log-likelihood should be returned.
#' @param gauss_kronrod_set Integer indicating how many nodes to use in Gauss-Kronrod quadrature. Values of 1-6 are supported and correspond to {11, 21, ..., 61} nodes, respectively.
#' @param gauss_leguerre_set Integer indicating how many nodes to use in Gauss-Leguerre quadrature. Values of 1-6 are supported and correspond to {10, 20, ..., 60} nodes, respectively.
#' @param n_cores Indicates the number of CPU threads to use.
#' @param d0 The starting diameter of tumors at onset.
#' 
#' @return Fitted Model.
#' \describe{
#'   \item{$onset}{The specified formula for the onset submodel.}
#'   \item{$growth}{The specified formula for the growth rate submodel.}
#'   \item{$sympt}{The specified formula for the symptomatic detection submodel.}
#'   \item{$sens}{The specified formula for the screening sensitivity submodel.}
#'   \item{$onset_dep_growth}{Indicates if the mean of the inverse growth rate distribution depends on the age at onset}
#'   \item{$par}{Numeric vector containing parameter values for all four submodels, either provided by start.par or assigned by the function.}
#'   \item{$fitted}{Indicates if the model parameters have been fitted to data. Set to TRUE if the optimization succeeds.}
#' }
#' @examples
#' 
#' @export
EstimateBCModel <- function(model, data, 
                            base_variables = list(
                              entry = "entry", exit = "exit", is_case = "case", 
                              mode = "mode", size = "size", scr_hist = "scr"), 
                            hessian = TRUE, 
                            gauss_kronrod_set = 6, gauss_leguerre_set = 4, 
                            n_cores = parallel::detectCores() - 1, d0 = 0.5) {
  bc_est_iter <<- 0
  cat("Estimating...\n")
  
  if (model$onset_dep_growth){
    IndL__ <- CalcIndividualLikelihoodODG
  } else {
    IndL__ <- CalcIndividualLikelihood
  }
  
  `%dopar%` <- foreach::`%dopar%`
  
  gkmat <- gknodes(gauss_kronrod_set)
  glmat <- glnodes(gauss_leguerre_set)
  
  v0 <- (pi / 6) * d0 ^ 3
  
  data <- MatchDataVariables(data, base_variables)
  
  CalcLogL <- function(par, data, d0) {
    bc_est_iter <<- bc_est_iter + 1
    new_par <- BuildParameters(par, model, data)
    
    log_lik <- foreach::foreach(i = itertools::isplitVector(1:nrow(data), chunks = n_cores), 
                                #.packages = "BCNatHist", 
                                .export = c("gkmat", "glmat", "IndL__", "v0"), 
                                .combine = sum) %dopar%
      vapply(i, function(i) {
        log(IndL__(data$case[i], data$mode[i], data$exit[i],
                   data$v[i], data$scr[[i]],
                   data$entry[i], data$e_scr[[i]],
                   new_par[i, ], gkmat$x, gkmat$w, glmat$x, glmat$w, d0 = d0, v0 = v0))
      }, FUN.VALUE = 0)
    cat("Iteration", bc_est_iter, ", time elapsed=", round(Sys.time() - start_time, 3), ": LogL=", log_lik, "    \r")
    ifelse(is.nan(log_lik), -1e16, log_lik)
  }
  
  start_time <- Sys.time()
  cl <- parallel::makePSOCKcluster(n_cores)
  doParallel::registerDoParallel(cl)
  est <- optim(model$par, function(pars) CalcLogL(pars, data, d0), method = 'BFGS',
               hessian = hessian, control = list(fnscale = -1))
  names(est$par) <- names(model$par)
  model$par <- est$par
  model$LogL <- est$value
  if (hessian) {
    model$hess <- est$hessian
  }
  parallel::stopCluster(cl)
  if (est$convergence == 0) {
    model$fitted <- TRUE
    model$info <- c(time = round(Sys.time() - start_time, 3), iterations = bc_est_iter)
    cat("\nEstimation successful! Maximum value found:", est$value, "\n")
  } else {
    cat("\nEstimation completed, but with message:\n", est$message, "\n")
  }
  return(model)
}
