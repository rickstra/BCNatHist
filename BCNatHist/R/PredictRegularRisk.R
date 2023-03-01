#' @title PredictRegularRisk
#' @description Predict BC risk, both if and if not attending the next screening.
#'
#' @param model Natural history model of class \code{BCModel} as defined by \code{DefineBCModel}.
#' @param data Dataframe containing the prediction data. One entry per individual.
#' @param base_variables List linking variable names of \code{data} to variables used in the prediction: 
#' \describe{
#'   \item{entry}{The age at the start of the prediction.}
#'   \item{exit}{The age at end of prediction follow-up.}
#'   \item{is_case}{BC case status (case = 1).}
#'   \item{mode}{Mode of detection (screen-detected = 1).}
#'   \item{size}{tumor size (mm).}
#'   \item{scr_hist}{Screening histories (each entry a list of ages screened).}
#' }
#' @param gauss_kronrod_set Integer indicating how many nodes to use in Gauss-Kronrod quadrature. Values of 1-6 are supported and correspond to {11, 21, ..., 61} nodes, respectively.
#' @param gauss_leguerre_set Integer indicating how many nodes to use in Gauss-Leguerre quadrature. Values of 1-6 are supported and correspond to {10, 20, ..., 60} nodes, respectively..
#' @param d0 The starting diameter of tumors at onset.
#' @param n_cores Indicates the number of CPU threads to use.
#'
#' @return \code{data.frame} with predictions of breast risk.
#' \describe{
#'   \item{$age}{Age at the start of the prediction interval.}
#'   \item{$end}{Age at the end of the prediction interval.}
#'   \item{$surv}{The model fitted probability of surviving until \code{age}.}
#'   \item{$risk_attend}{The predicted risk of diagnosis, if attending the next screening.}
#'   \item{$risk_noattend}{The predicted risk of diagnosis, if not attending the next screening.}
#' }
#' @examples
#' 
#' @export
PredictRegularRisk <- function(model, data, 
                                  base_variables = list(
                                    entry = "entry", exit = "exit", is_case = "case", 
                                    mode = "mode", size = "size", scr_hist = "scr"),
                                  gauss_kronrod_set = 6, gauss_leguerre_set = 4, 
                                  d0 = 0.5, n_cores = parallel::detectCores() - 1){
  
  `%dopar%` <- foreach::`%dopar%`
  
  gkmat <- gknodes(gauss_kronrod_set)
  glmat <- glnodes(gauss_leguerre_set)
  
  dv <- function(d) {
    pi / 6 * d ^ 3
  }
  v0 <- dv(d0)
  
  newpar <- BuildParameters(model$par, model, data)
  data <- MatchDataVariables(data, base_variables, include_prevalent_screen = FALSE)
  
  ix <- seq_along(data$entry)
  cl <- parallel::makePSOCKcluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  cat("Predicting...1/3    \r")
  # probability of surviving until start, i.e. left truncation
  q <- foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                        .combine = c) %dopar%
    vapply(i, function(j) {
      CalcCensCase(data$entry[j], data$e_scr[[j]], 
                   newpar[j, ], gkmat$x, gkmat$w, glmat$x, glmat$w)}, 
      FUN.VALUE = 0)
  cat("Predicting...2/3    \r")
  # risk of diagnosis if attending the next screening(s)
  overall_attend <- 1 - foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                                         .combine = c) %dopar%
    vapply(i, function(j) {
      CalcCensCase(data$exit[j], data$scr[[j]], 
                   newpar[j, ], gkmat$x, gkmat$w, glmat$x, glmat$w)}, 
      FUN.VALUE = 0) / q
  cat("Predicting...3/3    \n")
  # risk if *not* attending the next screening(s)
  overall_noattend <- 1 - foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                                           .combine = c) %dopar%
    vapply(i, function(j) {
      CalcCensCase(data$exit[j], data$e_scr[[j]], 
                   newpar[j, ], gkmat$x, gkmat$w, glmat$x, glmat$w)}, 
      FUN.VALUE = 0) / q
  parallel::stopCluster(cl)
  
  out <- data.frame(age = data$entry, end = data$exit, surv = q, 
                    risk_attend = overall_attend, risk_noattend = overall_noattend)
  cat("Prediction complete!\n")
  out
}