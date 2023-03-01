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

