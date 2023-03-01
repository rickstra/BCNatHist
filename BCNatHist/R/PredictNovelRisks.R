#' @title PredictNovelRisks
#' @description Predict novel BC risk, i.e. mode-specific risks and risk of onset.
#'
#' @param model Natural history model of class \code{BCModel} as defined by \code{DefineBCModel}.
#' @param data Dataframe containing the prediction data. One entry per individual.
#' @param regular_pred Optional dataframe containing the risk predictions from \code{PredictRegularRisk}. If not provided, they will be calculated first. 
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
#' @param followup Optional parameter--a positive number--setting the follow-up time for all participants. Supercedes the exit variable. 
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
PredictNovelRisks <- function(model, data, regular_pred = NULL,  
                              base_variables = list(
                                entry = "entry", exit = "exit", is_case = "case", 
                                mode = "mode", size = "size", scr_hist = "scr"),
                              gauss_kronrod_set = 6, gauss_leguerre_set = 4, followup = NULL,
                              d0 = 0.5, n_cores = parallel::detectCores() - 1){
  if(is.null(followup)) {
    followup <- data$exit
  } else {
    followup <- data$entry + followup
  }
  
  `%dopar%` <- foreach::`%dopar%`
  
  gkmat <- gknodes(gauss_kronrod_set)
  glmat <- glnodes(gauss_leguerre_set)
  
  dv <- function(d) {
    pi / 6 * d ^ 3
  }
  
  v_max <- dv(d_max)
  v0 <- dv(d0)
  
  newpar <- BCmodelParBuilder(model, data)
  data <- MatchDataVariables(data, base_variables, include_prevalent_screen = FALSE)
  if(is.null(followup)) {
    followup <- data$exit
  } else {
    followup <- data$entry + followup
  }
  
  ix <- seq_along(data$entry)
  parallel::cl <- makePSOCKcluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  if (is.null(regular_pred)) {
    cat("Predicting...1/7    \r")
    q <- foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                          .combine = c) %dopar% vapply(i, function(j) {
           CalcCensCase(data$entry[j],data$e_scr[[j]], newpar[j, ],
                        gkmat$x, gkmat$w, glmat$x, glmat$w)
           }, FUN.VALUE = 0)
    cat("Predicting...2/7    \r")
    overall <- 1 - foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                                    .combine = c) %dopar% vapply(i, function(j) {
                     CalcCensCase(followup[j], data$scr[[j]], newpar[j, ],
                                  gkmat$x, gkmat$w, glmat$x, glmat$w) 
                     }, FUN.VALUE = 0) / q
    cat("Predicting...3/7    \r")
    no_overall <- 1 - foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                                       .combine = c) %dopar% vapply(i, function(j) {
                        CalcCensCase(followup[j], data$e_scr[[j]], newpar[j, ],
                                     gkmat$x, gkmat$w, glmat$x, glmat$w) 
                        }, FUN.VALUE = 0) / q
  } else {
    q <- regular_pred$surv
    overall <- regular_pred$risk_attend
    no_overall <- regular_pred$risk_noattend
  }
  cat("Predicting...4/7    \r")
  onset0 <- pMVK(data$entry, newpar[, 1], newpar[, 2], newpar[, 3]) / q
  cat("Predicting...5/7    \r")
  onset <- pMVK(followup, newpar[, 1], newpar[, 2], newpar[, 3]) / q
  
  # screening now
  cat("Predicting...6/7    \r")
  sc0 <- foreach(i=isplitVector(ix, chunks = n_cores), .combine = c) %dopar%
    vapply(i, function(i) {
      CalcFutureScreeningSens(data$entry[i], data$e_scr[[i]][-1], newpar[i, ], 
                              gkmat$x, gkmat$w, glmat$x, glmat$w, d0, v0)
      }, FUN.VALUE = 0) / q
  # next planned screening
  cat("Predicting...7/7    \r")
  sc <- foreach(i=isplitVector(ix, chunks = n_cores), .combine = c) %dopar%
    vapply(i, function(i) {
      CalcFutureScreeningSens(data$scr[[i]][1], data$scr[[i]][-1], newpar[i, ], 
                              gkmat$x, gkmat$w, glmat$x, glmat$w, d0, v0)
      }, FUN.VALUE = 0) / q
  sym <- overall - sc
  #-----------------------------------------------------------------------------
  parallel::stopCluster(cl)
  out <- cbind(age = data$entry, end = followup, surv = q, 
               risk_overall = overall, risk_onset0 = ons0, 
               risk_onset = ons, screening_now = sc0, screening_next = sc, 
               risk_interval = sym)
  
  cat("Estimation successful!\n")
  out
}

