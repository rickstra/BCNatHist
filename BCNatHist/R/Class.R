#' @title DefineBCModel
#' @description Set up a natural history model.
#' @param onset  Specifies covariates to include in the onset submodel.
#' @param growth Specifies covariates to include in the inverse growth rate submodel.
#' @param sympt  Specifies covariates to include in the symptomatic detection submodel.
#' @param sens   Specifies covariates to include in the screening sensitivity submodel.
#' @param onset_dep_growth Indicate if the inverse growth rate should depend on the age at onset (Strandberg et al. 2022).
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
                          onset_dep_growth = FALSE, start_par = NULL) {
  var_names <- c(all.vars(onset), all.vars(growth), all.vars(sympt), all.vars(sens))
  dummy_df <- as.data.frame(matrix(1, 1, length(var_names)))
  names(dummy_df) <- var_names
  
  o_names <- colnames(model.matrix(onset, dummy_df))[-1]
  g_names <- colnames(model.matrix(growth, dummy_df))[-1]
  sym_names <- colnames(model.matrix(sympt, dummy_df))[-1]
  sen_names <- colnames(model.matrix(sens, dummy_df))[-1]
  
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

#===============================================================================

#' @title EstimateBCModel
#' @description Fit a BCModel to cohort data.
#'
#' @param model Natural history model of class \code{BCModel} as defined by \code{DefineBCModel}.
#' @param data Dataframe containing the cohort data. One entry per individual.
#' @param hessian Indicates if the Hessian matrix of the maximum log-likelihood should be returned.
#' @param n_cores Indicates the number of CPU threads to use.
#' @param d0 The starting diameter of tumors at onset.
#' @param entry The variable name in \code{data} corresponding to the age at study entry.
#' @param exit  The variable name in \code{data} corresponding to the age at end of follow-up. For BC cases, this is also the age at diagnosis. 
#' @param is_case The variable name in \code{data} indicating BC case status (case = 1).
#' @param mode The variable name in \code{data} indicating the mode of detection (screen-detected = 1).
#' @param size The variable name in \code{data} corresponding to the tumor size (mm).
#' @param scr_hist The variable name in \code{data} corresponding to the screening histories (each entry a list of ages screened).
#' @param gauss_kronrod_set Integer indicating how many nodes to use in Gauss-Kronrod quadrature. Values of 1-6 are supported and correspond to {11, 21, ..., 61} nodes, respectively.
#' @param gauss_leguerre_set Integer indicating how many nodes to use in Gauss-Leguerre quadrature. Values of 1-6 are supported and correspond to {10, 20, ..., 60} nodes, respectively..
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

EstimateBCModel <- function(model, data, hessian = TRUE, 
                            n_cores = parallel::detectCores() - 1, d0 = 0.5, 
                            entry = "entry", exit = "exit", is_case = "case", 
                            mode = "mode", size = "size", scr_hist = "scr", 
                            gauss_kronrod_set = 6, gauss_leguerre_set = 4) {
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
  
  data$entry <- data[, entry]
  data$exit <- data[, exit]
  data$case <- data[, is_case]
  data$mode <- data[, mode]
  data$v <- (pi / 6) * data[, size] ^ 3
  data$scr <- data[, scr_hist]
  
  data$e_scr <- lapply(1:nrow(data), function(i) {
      data$scr[[i]][data$scr[[i]] < data$entry[i]]
    })
  start_time <- Sys.time()
  
  BuildParams <- function(par, data, model) {
    
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
    # [0] A
    # [1] B
    # [2] delta
    # [3] phi
    # [4] mu
    # [5] eta
    # [6] beta_1
    # [7] beta_0
    # [8] odg (optional)
  }
  
  CalcLogL <- function(par, data, d0) {
    bc_est_iter <<- bc_est_iter + 1
    new_par <- BuildParams(par, data, model)
    
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

#' @title SummarizeBCEstimates
#' @description Description
#'
#' @param model Natural history model of class \code{BCModel} as defined by \code{DefineBCModel}.
#' @param alpha Par.
#' @param trans Par.
#' 
#' @return Fitted Model.
#' @export
SummarizeBCEstimates <- function(model, alpha = 0.05, trans = TRUE) {
  with(model, {
    if (!is.null(model$hess)){
    cov <- solve(-hess)
    }
    se <- sqrt(diag(cov))
    q <- qnorm(1 - alpha/2)
    p <- round(2 * (1 - pnorm(abs(as.numeric(par) / se))), 3)
    sig <- paste0(ifelse(p <= 0.1, '*', ''), ifelse(p <= 0.05, '*', ''), 
                  ifelse(p <= 0.01, '*', ''), ifelse(p <= 0.005, '*', '')) 
    sig[names(par) %in% c('A', 'B', 'phi', 'beta1', 'delta0', 'mu0', 'eta0', 'beta0')] <- 'NA'
    p[names(par) %in% c('A', 'B', 'phi', 'beta1', 'delta0', 'mu0', 'eta0', 'beta0')] <- NA
    out <- data.frame(parameter = names(par), estimate = par, lowerCI = par - q * se, upperCI = par + q * se, p = p, sig)
    out[1, 2:4] <- -exp(out[1, 2:4])
    out[2:3, 2:4] <- exp(out[2:3, 2:4])
    if (trans) {
      ix <- which(out$par == 'beta0')
      out[5:(ix - 1), 2:4] <- exp(out[5:(ix - 1), 2:4])
      if (onset_dep_growth){
        out[nrow(out), 2:4] <- exp(out[nrow(out), 2:4])
      }
    }
    rownames(out) <- NULL
    out[, 2:4] <- round(out[, 2:4], 6)
    return(out)
  })
}

#' @title summary.BCModel
#' @description Description
#'
#' @param model Natural history model of class \code{BCModel} as defined by \code{DefineBCModel}.
#' @return Fitted Model.
#' @export
summary.BCModel <- function(model, ...) {
  cat("Call:\n")
  cat("Onset:  ")
  print(model$onset)
  cat("Growth: ")
  print(model$growth)
  cat("Sympt:  ")
  print(model$sympt)
  cat("Sens:   ")
  print(model$sens)
  cat("Onset-dependent growth: ")
  print(model$onset_dep_growth)
  cat("#####################################################\n")
  if (model$fitted) {
    print(SummarizeBCEstimates(model, ...))
  }
  invisible(SummarizeBCEstimates(model, ...))
}

