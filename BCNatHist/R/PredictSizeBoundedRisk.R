#' @title PredictSizeBoundedRisk
#' @description Predict BC risk with an upper limit on tumor size, both if and if not attending the next screening.
#'
#' @param model Natural history model of class \code{BCModel} as defined by \code{DefineBCModel}.
#' @param data Dataframe containing the cohort data. One entry per individual.
#' @param d_max The upper limit for the tumor diameter.
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
#'   \item{$bounded_attend}{The predicted risk of diagnosis with a tumor less than \code{d_max}, if attending the next screening.}
#'   \item{$bounded_noattend}{The predicted risk of diagnosis with a tumor less than \code{d_max}, if not attending the next screening.}
#' }
#' age = data$entry, end = data$exit, surv = q, 
# risk_attend = overall, risk_noattend = no_overall,
# bounded_attend = d10_int1 + d10_scr + d10_int2scr, 
# bounded_noattend = d10_int1 + d10_int2noscr)
#' @examples
#' 
#' @export
PredictSizeBoundedRisk <- function(model, data, d_max = 10, regular_pred = NULL,  
                                  base_variables = list(
                                    entry = "entry", exit = "exit", is_case = "case", 
                                    mode = "mode", size = "size", scr_hist = "scr"),
                                  gauss_kronrod_set = 6, gauss_leguerre_set = 4, followup = NULL,
                                  d0 = 0.5, n_cores = parallel::detectCores() - 1){
  
  `%dopar%` <- foreach::`%dopar%`
  
  gkmat <- gknodes(gauss_kronrod_set)
  glmat <- glnodes(gauss_leguerre_set)
  
  dv <- function(d) {
    pi / 6 * d ^ 3
  }
  
  v_max <- dv(d_max)
  v0 <- dv(d0)
  
  newpar <- BuildParameters(model$par, model, data)
  data <- MatchDataVariables(data, base_variables, include_prevalent_screen = FALSE)
  if(is.null(followup)) {
    followup <- data$exit
  } else {
    followup <- data$entry + followup
  }
  
  
  ix <- seq_along(data$entry)
  cl <- parallel::makePSOCKcluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  if (is.null(regular_pred)) {
    cat("Predicting...1/7    \r")
    q <- foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                          .combine = c) %dopar%
      vapply(i, function(j) {
        CalcCensCase(data$entry[j],data$e_scr[[j]], newpar[j, ],
                     gkmat$x, gkmat$w, glmat$x, glmat$w)}, FUN.VALUE = 0)
    cat("Predicting...2/7    \r")
    overall <- 1 - foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                                    .combine = c) %dopar% {
                                      vapply(i, function(j) {
                                        CalcCensCase(followup[j], data$scr[[j]], newpar[j, ],
                                                     gkmat$x, gkmat$w, glmat$x, glmat$w) 
                                      }, FUN.VALUE = 0)} / q
    cat("Predicting...3/7    \r")
    no_overall <- 1 - foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                                       .combine = c) %dopar% {
                                         vapply(i, function(j) {
                                           CalcCensCase(followup[j], data$e_scr[[j]], newpar[j, ],
                                                        gkmat$x, gkmat$w, glmat$x, glmat$w) 
                                         }, FUN.VALUE = 0)} / q
  } else {
    q <- regular_pred$surv
    overall <- regular_pred$risk_attend
    no_overall <- regular_pred$risk_noattend
  }
  
  cat("Predicting...4/7    \r")
  dmax_scr <- foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                              .combine = c) %dopar% {
                                dj.grid <- expand.grid(d=gkmat$x * (d_max - d0) / 2 + (d_max + d0) / 2,
                                                       j=i)
                                temp <- mapply(function(d, j) {pi * d ^ 2 / 2 * 
                                    CalcScreenCase(data$scr[[j]][1], dv(d), data$scr[[j]], newpar[j, ], gkmat$x, gkmat$w)},
                                    dj.grid$d, dj.grid$j)
                                matrix(temp, , nrow(gkmat), byrow = T) %*% (gkmat$w * (d_max - d0) / 2)
                                #temp
                              } /q
  cat("Predicting...5/7    \r")
  dmax_int1 <- foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                               .combine = c) %dopar% {
                                 vapply(i, function(j) {
                                   ud_grid <- expand.grid(u = (gkmat$x * (data$scr[[j]][1] - data$entry[j]) + 
                                                                 data$scr[[j]][1] + data$entry[j]) / 2,
                                                          d = (gkmat$x * (d_max - d0) + d_max + d0) / 2)
                                   wts <- (gkmat$w * (d_max - d0) / 2) %*% 
                                     t((data$scr[[j]][1] - data$entry[j]) / 2 * gkmat$w)
                                   temp <- mapply(function(u, d) {
                                     CalcSymptCase(u, dv(d), data$e_scr[[j]], newpar[j, ], 
                                                   gkmat$x, gkmat$w) * pi * d ^ 2 / 2 
                                   }, ud_grid$u, ud_grid$d)
                                   sum(matrix(temp, nrow(wts), , T) * wts)
                                 }, FUN.VALUE = 0)} / q
  cat("Predicting...6/7    \r")
  dmax_int2scr <- foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores), 
                                  .combine = c) %dopar% {
                                    vapply(i, function(j) {
                                      ud_grid <- expand.grid(u = (gkmat$x * (followup[j] - data$scr[[j]][1]) + 
                                                                    data$scr[[j]][1] + followup[j]) / 2,
                                                             d = (gkmat$x * (d_max - d0) + d_max + d0) / 2)
                                      wts <- (gkmat$w * (d_max - d0) / 2) %*% 
                                        t((followup[j] - data$scr[[j]][1]) / 2 * gkmat$w)
                                      temp <- mapply(function(u, d) {
                                        CalcSymptCase(u, dv(d), data$scr[[j]], newpar[j, ], 
                                                      gkmat$x, gkmat$w) * pi * d ^ 2 / 2 
                                      }, ud_grid$u, ud_grid$d)
                                      sum(matrix(temp, nrow(wts), , T) * wts)
                                    }, FUN.VALUE = 0)} / q
  cat("Predicting...7/7    \n")
  dmax_int2noscr <- foreach::foreach(i=itertools::isplitVector(ix, chunks = n_cores),
                                    .combine = c) %dopar% {
                                      vapply(i, function(j) {
                                        ud_grid <- expand.grid(u = (gkmat$x * (followup[j] - data$scr[[j]][1]) + 
                                                                      data$scr[[j]][1] + followup[j]) / 2,
                                                               d = (gkmat$x * (d_max - d0) + d_max + d0) / 2)
                                        wts <- (gkmat$w * (d_max - d0) / 2) %*% 
                                          t((followup[j] - data$scr[[j]][1]) / 2 * gkmat$w)
                                        temp <- mapply(function(u, d) {
                                          CalcSymptCase(u, dv(d), data$e_scr[[j]], newpar[j, ], gkmat$x, 
                                                        gkmat$w) * pi * d ^ 2 / 2 
                                        }, ud_grid$u, ud_grid$d)
                                        sum(matrix(temp, nrow(wts), , T) * wts)
                                      }, FUN.VALUE = 0)}  / q
  parallel::stopCluster(cl)
  cat("Prediction successful!\n")
  data.frame(age = data$entry, end = data$exit, surv = q, 
             risk_attend = overall, risk_noattend = no_overall,
             bounded_attend = dmax_int1 + dmax_scr + dmax_int2scr, 
             bounded_noattend = dmax_int1 + dmax_int2noscr)
}
