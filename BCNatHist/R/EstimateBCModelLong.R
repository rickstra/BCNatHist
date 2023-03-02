#' @title EstimateBCModelLong
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
EstimateBCModelLong <- function(model, data, scr_data, 
                                link = NULL, e_link = NULL, 
                                age_mammogram = "age_mammogram", id = "id",
                                base_variables = list(
                                  entry = "entry", exit = "exit", is_case = "case", 
                                  mode = "mode", size = "size", scr_hist = "scr"), 
                                hessian = TRUE, gauss_kronrod_set = 6, gauss_leguerre_set = 4, 
                                n_cores = parallel::detectCores() - 1, d0 = 0.5) {
  bc_est_iter__ <<- 0
  cat("Estimating...\n")
  
  IndL__ <- CalcIndividualLikelihoodLong
  
  `%dopar%` <- foreach::`%dopar%`
  
  gkmat <- gknodes(gauss_kronrod_set)
  glmat <- glnodes(gauss_leguerre_set)
  
  v0 <- (pi / 6) * d0 ^ 3
  
  if (is.null(link) | is.null(e_link)) {
    link <- lapply(1:nrow(data), function(i) which(scr_data[, id] == data[i, id]))
    e_link <- lapply(1:nrow(data), function(i) which(scr_data[, id] == data[i, id] & 
                  scr_data[, age_mammogram] < data[i, ifelse(is.null(base_variables$entry), "entry", base_variables$entry)]))
  } else if (sum(sapply(link, function(lst) length(unique(p4data$studieid[lst]))) != 1) != 0) {
    print("Link mismatch between base data and screening data")
    Wah_wah()
  }
  
  data <- MatchDataVariables(data, base_variables)
  
  g <- function(par, data, scr_data, link, e_link, d0) {
    bc_est_iter__ <<- bc_est_iter__ + 1
    newpar <- ParBuilder(par, data, scr_data, model)
    
    logL <- foreach::foreach(i=itertools::isplitVector(1:nrow(data), chunks = n.cores), 
                    .combine = sum) %dopar%
      vapply(i, function(i) {
        log(IndL__(data$case[i], data$mode[i], data$exit[i],
                   data$v[i], scr_data$age_mammogram[link[[i]]],
                   data$entry[i], scr_data$age_mammogram[e_link[[i]]],
                   newpar[[1]][i, ], newpar[[2]][link[[i]]], newpar[[2]][e_link[[i]]],
                   gkmat$x, gkmat$w, glmat$x, glmat$w, d0 = d0, v0 = v0))
      }, FUN.VALUE = 0)
    cat("Iteration", bc_est_iter__, ", time elapsed=", 
        round(Sys.time() - start.time, 3), ": LogL=", logL, '    \r')
    ifelse(is.nan(logL)|is.infinite(logL), -1e16, logL)
  }
  
  cl <- parallel::makePSOCKcluster(n_cores)
  doParallel::registerDoParallel(cl)
  start.time <- Sys.time()
  est <- optim(model$par, function(pars) g(pars, data, scr_data, link, e_link, d0), 
               method = 'BFGS', hessian = hessian, control = list(fnscale = -1))
  names(est$par) <- names(model$par)
  model$par <- est$par
  model$LogL <- est$value
  if (hessian) {
    model$hess <- est$hessian
  }
  parallel::stopCluster(cl)
  if (est$convergence == 0) {
    model$fitted = TRUE
    model$info <- c(time = round(Sys.time() - start.time, 3), iterations = bc_est_iter__)
    cat("\nEstimation successful! Maximum value found:", est$value, '\n')
  } else {
    cat("\nEstimation completed, but with message:\n", est$message, '\n')
  }
  model
}


