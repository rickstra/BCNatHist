
###################################################
# DefineBCmodel: Set up a natural history model
###################################################
# Arguments
# onset:  specifies covariates to include in the onset submodel.
# growth: specifies covariates to include in the inverse growth rate submodel.
# sympt:  specifies covariates to include in the symptomatic detection submodel.
# sens:   specifies covariates to include in the screening sensitivity submodel.
# onset.dep.growth: indicate if the inverse growth rate should depend on the 
#                   age at onset (Study 2).
# start_par: optionally specify (starting) parameter values.
#####
# Value
# A list-type object of class "BCmodel".
# onset:    the covariates of the onset submodel.
# growth:   the covariates of the inverse growth rate submodel.
# sympt:    the covariates of the symptomatic detection submodel.
# sens:     the covariates of the screening sensitivity submodel.
# OnsetDependentGrowth: indicates if the inverse growth rate depends on the 
#                       age at onset (Study 2).
# par:      numeric vector containing parameter values for all four submodels, 
#           either provided by start_par or assigned by the function.
# fitted:   indicates if the model has been fitted to data. Set to FALSE.
DefineBCmodel <- function(onset = ~ 1, growth = ~ 1, sympt = ~ 1, sens = ~ 1, 
                          onset.dep.growth = FALSE, free.size = FALSE, start_par = NULL) {
  vnames <- c(all.vars(onset), all.vars(growth), all.vars(sympt), all.vars(sens))
  xd <- as.data.frame(matrix(1, 1, length(vnames)))
  names(xd) <- vnames
  
  o_names <- colnames(model.matrix(onset, xd))[-1]
  g_names <- colnames(model.matrix(growth, xd))[-1]
  sym_names <- colnames(model.matrix(sympt, xd))[-1]
  sen_names <- colnames(model.matrix(sens, xd))[-1]
  
  parnames <- c('A', 'B', 'phi', 'beta1', 'delta0', o_names, 'mu0', g_names, 'eta0', sym_names, 'beta0', sen_names)
  
  if(is.null(start_par)){
    start_par <- c(-2.28, -7.09, -0.49, 0.47,
                   -3.08, sample(c(-0.01, 0.01), length(o_names), T),
                   0.07, sample(c(-0.01, 0.01), length(g_names), T),
                   -8.85, sample(c(-0.01, 0.01), length(sym_names), T),
                   -5.37, sample(c(-0.01, 0.01), length(sen_names), T))
  }
  if (onset.dep.growth) {
    parnames <- c(parnames, 'odg')
    start_par <- c(start_par, 0.01)
  } else if (free.size) {
    parnames <- c(parnames, 'z')
    start_par <- c(start_par, 1)
  }
  names(start_par) <- parnames
  model <- list(onset = onset, growth = growth, sympt = sympt, sens = sens, 
                "OnsetDependentGrowth" = onset.dep.growth, "FreeSize" = free.size,
                par = start_par, fitted = FALSE)
  class(model) <- "BCmodel"
  model
}
#####
DefineBCmodel2 <- function(data = NULL, scr_data = NULL, onset = ~ 1, growth = ~ 1, sympt = ~ 1, sens = ~ 1, 
                          onset.dep.growth = FALSE, start_par = NULL) {
  if (is.null(data)) {
    vnames <- c(all.vars(onset), all.vars(growth), all.vars(sympt), all.vars(sens))
    xd <- as.data.frame(matrix(1, 1, length(vnames)))
    names(xd) <- vnames
    
    
    o_names <- colnames(model.matrix(onset, xd))[-1]
    g_names <- colnames(model.matrix(growth, xd))[-1]
    sym_names <- colnames(model.matrix(sympt, xd))[-1]
    sen_names <- colnames(model.matrix(sens, xd))[-1]
  }
  
  parnames <- c('A', 'B', 'phi', 'beta1', 'delta0', o_names, 'mu0', g_names, 'eta0', sym_names, 'beta0', sen_names)
  
  # model.matrix(~ ., data=testFrame, 
  #              contrasts.arg = lapply(testFrame[, sapply(testFrame, is.factor)], contrasts, contrasts=FALSE))
  
  
  if(is.null(start_par)){
    start_par <- c(-2.28, -7.09, -0.49, 0.47,
                   -3.08, sample(c(-0.01, 0.01), length(o_names), T),
                   0.07, sample(c(-0.01, 0.01), length(g_names), T),
                   -8.85, sample(c(-0.01, 0.01), length(sym_names), T),
                   -5.37, sample(c(-0.01, 0.01), length(sen_names), T))
  }
  if (onset.dep.growth) {
    parnames <- c(parnames, 'odg')
    start_par <- c(start_par, 0.01)
  }
  names(start_par) <- parnames
  model <- list(onset = onset, growth = growth, sympt = sympt, sens = sens, 
                "OnsetDependentGrowth" = onset.dep.growth, 
                par = start_par, fitted = FALSE)
  class(model) <- "BCmodel"
  model
}

EstimateBCmodel <- function(model, data, logsize = FALSE, hessian = TRUE, 
                            n.cores = parallel::detectCores() - 1, d0 = 0.5) {
  library(foreach)
  library(doParallel)
  library(itertools)
  bc_est_iter__ <<- 0
  cat("Estimating...\n")
  
  if (model$OnsetDependentGrowth){
    IndL__ <- BC::IndL_odg
  } else if (model$FreeSize) {
    IndL__ <- BC::IndL2
  } else {
    IndL__ <- BC::IndL
  }
  
  ParBuilder <- function(par, data, model) {
    
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
    
    parnames <- c('A', 'B', 'delta', 'phi', 'mu', 'eta', 'beta1', 'beta2')
    # Onset-dependent growth?
    if(model$OnsetDependentGrowth) {
      parnames <- c(parnames, "odg")
    # } else if (model$FreeSize) {
    #   parnames <- c(parnames, "z")
    } else {
      par <- NULL
    }
    
    par.mat <- cbind(A, B, exp(o_pred), phi, exp(g_pred), 
                 exp(sym_pred), beta, scr_pred, par)
    colnames(par.mat) <- parnames
    par.mat
    # [0] A
    # [1] B
    # [2] delta
    # [3] phi
    # [4] mu
    # [5] eta
    # [6] beta_1
    # [7] beta_0
    # [8] odg or z (optional)
  }
  g <- function(par, data, d0) {
    bc_est_iter__ <<- bc_est_iter__ + 1
    newpar <- ParBuilder(par, data, model)
    
    ifelse(model$FreeSize, par[length(par)], "")
    
    logL <- foreach(i=isplitVector(1:nrow(data), chunks = n.cores), .packages = "BC", .export = c("gkmat61", "glmat40", "IndL__"), .combine = sum) %dopar%
      sapply(i, function(i) {
        log(IndL__(data$case[i], data$mode[i], data$exit[i],
                     data$v[i], data$scr[[i]],
                     data$entry[i], data$e.scr[[i]],
                     newpar[i, ], gkmat61$x, gkmat61$w, glmat40$x, glmat40$w, d0 = d0))
      })
    cat("Iteration", bc_est_iter__, ", time elapsed=", round(Sys.time() - start.time, 3), ": LogL=", logL, ifelse(model$FreeSize, round(par[length(par)], 2), ""), '    \r')
    ifelse(is.nan(logL), -1e16, logL)
  }
  
  cl <- makePSOCKcluster(n.cores)
  registerDoParallel(cl)
  start.time <- Sys.time()
  est <- optim(model$par, function(pars) g(pars, data, d0), method = 'BFGS',
               hessian = hessian, control = list(fnscale = -1))
  names(est$par) <- names(model$par)
  model$par <- est$par
  model$LogL <- est$value
  if (hessian) {
    model$hess <- est$hessian
  }
  stopCluster(cl)
  if (est$convergence == 0) {
    model$fitted = TRUE
    model$info <- c(time = round(Sys.time() - start.time, 3), iterations = bc_est_iter__)
    cat("\nEstimation successful! Maximum value found:", est$value, '\n')
  } else {
    cat("\nEstimation completed, but with message:\n", est$message, '\n')
  }
  #rm(bc_est_iter__)
  model
  
  # ParBuilder(model$par, data, model$onset, model$growth, model$sympt, model$sens)
}

BCmodelCI <- function(model, alpha = 0.05, trans = FALSE) {
  with(model, {
    if (!is.null(model$hess)){
    cov <- solve(-hess)
    }
    se <- sqrt(diag(cov))
    q <- qnorm(1 - alpha/2)
    p <- round(2*(1-pnorm(abs(as.numeric(par)/se))), 3)
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
      if (OnsetDependentGrowth){
        out[nrow(out), 2:4] <- exp(out[nrow(out), 2:4])
      }
    }
    rownames(out) <- NULL
    out[, 2:4] <- round(out[, 2:4], 6)
    return(out)
  })
}

summary.BCmodel <- function(x, ...) {
  cat("Call:\n")
  cat("Onset:  ")
  print(x$onset)
  cat("Growth: ")
  print(x$growth)
  cat("Sympt:  ")
  print(x$sympt)
  cat("Sens:   ")
  print(x$sens)
  cat("Onset-dependent growth: ")
  print(x$OnsetDependentGrowth)
  cat("#####################################################\n")
  if (x$fitted) {
    print(BCmodelCI(x, ...))
  }
  invisible(BCmodelCI(x, ...))
}

###########################################
## Other stuff
###########################################

BCmodelParBuilder <- function(model, data) {
  
  # linear predictors
  o_pred <- model.matrix(model$onset, data)
  g_pred <- model.matrix(model$growth, data)
  sym_pred <- model.matrix(model$sympt, data)
  scr_pred <- model.matrix(model$sens, data)
  
  # Core
  A <- -exp(model$par[1])
  B <- exp(model$par[2])
  phi <- exp(model$par[3])
  beta <- model$par[4]
  par <- model$par[-(1:4)]
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
  
  # par <- cbind(exp(sym_pred), exp(g_pred), phi, scr_pred, beta, A, B, exp(o_pred))
  par <- cbind(A, B, delta = exp(o_pred), phi = phi, mu = exp(g_pred), eta = exp(sym_pred), beta,  beta2 = scr_pred)
  colnames(par) <- c('A', 'B', 'delta', 'phi', 'mu', 'eta', 'beta1',  'beta')
  par
  # [0] A
  # [1] B
  # [2] delta
  # [3] phi
  # [4] mu
  # [5] eta
  # [6] beta_1
  # [7] beta_0
}

BCmodelRisk <- function(model, data, x, n.cores = parallel::detectCores() - 1) {
  library(foreach)
  library(doParallel)
  library(itertools)
  cat("Estimating...\n")
  newpar <- BCmodelParBuilder(model, data)
  ix <- which(data$case == 0)
  cl <- makePSOCKcluster(n.cores)
  registerDoParallel(cl)
      q <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                   .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
        sapply(i, function(i) {BC::f_CensCase(data$exit[i],data$scr[[i]], newpar[i, ],
                          gkmat61$x, gkmat61$w, glmat40$x, glmat40$w)})
      ons0 <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                     .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
        sapply(i, function(i) {BC::f_pMVK(data$exit[i], newpar[i, 1], newpar[i, 2], 
                                          newpar[i, 3])}) / q
      ons <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                     .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
        sapply(i, function(i) {BC::f_pMVK(data$exit[i] + x, newpar[i, 1], newpar[i, 2], 
                newpar[i, 3])}) / q
      sym <- 1 - foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                         .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
        sapply(i, function(i) {BC::f_CensCase(data$exit[i] + x, data$scr[[i]], newpar[i, ], 
                                gkmat61$x, gkmat61$w, glmat40$x, glmat40$w)}) / q
      sc0 <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                     .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
        sapply(i, function(i) {BC::f_RiskScreen(0, data$exit[i], data$scr[[i]], newpar[i, ],
                           gkmat61$x, gkmat61$w, glmat40$x, glmat40$w)}) / q
      sc <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                    .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
        sapply(i, function(i) {BC::f_RiskScreen(x, data$exit[i], data$scr[[i]], newpar[i, ],
                           gkmat61$x, gkmat61$w, glmat40$x, glmat40$w)}) / q
      eith <- 1 - foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                          .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
        sapply(i, function(i) {BC::f_CensCase(data$exit[i] + x, 
                    c(data$scr[[i]], data$exit[i] + x), newpar[i, ], gkmat61$x, 
                    gkmat61$w, glmat40$x, glmat40$w)}) / q
  stopCluster(cl)
  out <- cbind(data$studieid[ix], data$exit[ix], q, ons0, ons, sc0, sc, sym, eith)
  colnames(out) <- c('studieid', 'age', 'surv', paste0('ons', c(0, x)),
    paste0('scr', c(0, x)), paste0(c('sym', 'risk'), x)
    )
  cat("Estimation successful!\n")
  out
  
}

BCmodelRiskPred <- function(model, data, d_max = 10, n.cores = parallel::detectCores() - 1) {
  library(foreach)
  library(doParallel)
  library(itertools)
  
  ###
  d10_int1 <- 0
  d10_int2noscr <- 0
  d10_scr <- 0
  d10_int2scr <- 0
  ###
  
  RiskScreen__ <- ifelse(model$OnsetDependentGrowth, BC::f_RiskScreen_odg, BC::f_RiskScreen)
  CensCase__ <- ifelse(model$OnsetDependentGrowth, BC::f_CensCase_odg, BC::f_CensCase)
  ScreenCase__ <- ifelse(model$OnsetDependentGrowth, BC::f_ScreenCase_odg, BC::f_ScreenCase)
  SymptCase__ <- ifelse(model$OnsetDependentGrowth, BC::f_SymptCase_odg, BC::f_SymptCase)
  
  cat("Predicting...\n")
  newpar <- BCmodelParBuilder(model, data)
  ix <- seq_along(data$entry)
  cl <- makePSOCKcluster(n.cores)
  registerDoParallel(cl)
  q <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
               .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
    sapply(i, function(i) {CensCase__(data$entry[i],data$e.scr[[i]], newpar[i, ],
                                          gkmat61$x, gkmat61$w, glmat40$x, glmat40$w)})
  sc <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
    sapply(i, function(i) {
      sum(sapply(data$scr[[i]][-which(data$scr[[i]] %in% data$e.scr[[i]])], function(k) {
        temp.scr <- data$scr[[i]][which(data$scr[[i]] <= k)]
        RiskScreen__(k, temp.scr, 
                         newpar[i, ], gkmat61$x, 
                         gkmat61$w, glmat40$x, glmat40$w)
      }))
    }) / q
  screensize <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
          .packages = "BC", .export = c("gkmat61", "glmat40", "d0", "dv")) %dopar% {
    jd.grid <- expand.grid(d=glmat40$x+d0, j=i)
    mapply(function(j, d) {pi * d^2 / 2 * ScreenCase__(data$scr[[j]][1], dv(d), data$scr[[j]], newpar[j, ], gkmat61$x, gkmat61$w) / q[j]},
           jd.grid$j, jd.grid$d)}
    #matrix(temp, length(i), , F)

  intsize1 <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                       .packages = "BC", .export = c("gkmat61", "glmat40", "d0", "dv")) %dopar% {
                         sapply(i, function(j) {
                           ud.grid <- expand.grid(u=gkmat61$x * (data$scr[[j]][1]-data$entry[j])/2 + (data$scr[[j]][1]+data$entry[j])/2,
                                                  d=glmat40$x + d0)
                           wts <- rep(1, length(glmat40$x)) %*% t((data$scr[[j]][1]-data$entry[j])/2 * gkmat61$w)
                           temp <- mapply(function(u, d) {SymptCase__(u, dv(d), data$e.scr[[j]], newpar[j, ], gkmat61$x, gkmat61$w) * pi * d^2 / 2  / q[j] },
                                          ud.grid$u, ud.grid$d)
                           rowSums(matrix(temp, nrow(wts), , T) * wts, na.rm = T)
                         })}
  intsize2 <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                       .packages = "BC", .export = c("gkmat61", "glmat40", "d0", "dv")) %dopar% {
                         sapply(i, function(j) {
                           ud.grid <- expand.grid(u=gkmat61$x * (data$exit[j]-data$scr[[j]][1])/2 + (data$exit[j]+data$scr[[j]][1])/2,
                                                  d=glmat40$x + d0)
                           wts <- rep(1, length(glmat40$x)) %*% t((data$exit[j]-data$scr[[j]][1])/2 * gkmat61$w)
                           temp <- mapply(function(u, d) {SymptCase__(u, dv(d), data$scr[[j]], newpar[j, ], gkmat61$x, gkmat61$w) * pi * d^2 / 2  / q[j] },
                                          ud.grid$u, ud.grid$d)
                           rowSums(matrix(temp, nrow(wts), , T) * wts, na.rm = T)
                         })}
  ###################################################

  # d10_int1 <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #     .packages = "BC", .export = c("gkmat61", "glmat40", "d0", "dv")) %dopar% {
  #       sapply(i, function(j) {
  #         ud.grid <- expand.grid(u=gkmat61$x * (data$scr[[j]][1]-data$entry[j])/2 + (data$scr[[j]][1]+data$entry[j])/2,
  #                                d=gkmat61$x * (d_max - d0)/2 + (d_max + d0)/2)
  #         wts <- (gkmat61$w * (d_max-d0)/2) %*% t((data$scr[[j]][1]-data$entry[j])/2 * gkmat61$w)
  #         temp <- mapply(function(u, d) {SymptCase__(u, dv(d), data$e.scr[[j]], newpar[j, ], gkmat61$x, gkmat61$w) * pi * d^2 / 2 },
  #                ud.grid$u, ud.grid$d)
  #         sum(matrix(temp, nrow(wts), , T) * wts)
  #       })}  / q
  # d10_int2scr <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                      .packages = "BC", .export = c("gkmat61", "glmat40", "d0", "dv")) %dopar% {
  #                        sapply(i, function(j) {
  #                          ud.grid <- expand.grid(u=gkmat61$x * (data$exit[j]-data$scr[[j]][1])/2 + (data$scr[[j]][1]+data$exit[j])/2,
  #                                                 d=gkmat61$x * (d_max - d0)/2 + (d_max + d0)/2)
  #                          wts <- (gkmat61$w * (d_max-d0)/2) %*% t((data$exit[j]-data$scr[[j]][1])/2 * gkmat61$w)
  #                          temp <- mapply(function(u, d) {SymptCase__(u, dv(d), data$scr[[j]], newpar[j, ], gkmat61$x, gkmat61$w) * pi * d^2 / 2 },
  #                                         ud.grid$u, ud.grid$d)
  #                          sum(matrix(temp, nrow(wts), , T) * wts)
  #                        })}  / q
  # d10_int2noscr <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                         .packages = "BC", .export = c("gkmat61", "glmat40", "d0", "dv")) %dopar% {
  #                           sapply(i, function(j) {
  #                             ud.grid <- expand.grid(u=gkmat61$x * (data$exit[j]-data$scr[[j]][1])/2 + (data$scr[[j]][1]+data$exit[j])/2,
  #                                                    d=gkmat61$x * (d_max - d0)/2 + (d_max + d0)/2)
  #                             wts <- (gkmat61$w * (d_max-d0)/2) %*% t((data$exit[j]-data$scr[[j]][1])/2 * gkmat61$w)
  #                             temp <- mapply(function(u, d) {SymptCase__(u, dv(d), data$e.scr[[j]], newpar[j, ], gkmat61$x, gkmat61$w) * pi * d^2 / 2 },
  #                                            ud.grid$u, ud.grid$d)
  #                             sum(matrix(temp, nrow(wts), , T) * wts)
  #                           })}  / q
  # d10_scr <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                    .packages = "BC", .export = c("gkmat61", "glmat40", "d0", "dv")) %dopar% {
  #   dj.grid <- expand.grid(d=gkmat61$x * (10 - d0)/2 + (10 + d0)/2,
  #                          j=i)
  #   temp <- mapply(function(d, j) {pi * d^2 / 2 * ScreenCase__(data$scr[[j]][1], dv(d), data$scr[[j]], newpar[j, ], gkmat61$x, gkmat61$w)},
  #               dj.grid$d, dj.grid$j)
  #   matrix(temp, , nrow(gkmat61), byrow = T) %*% (gkmat61$w * (10 - d0)/2)
  #   #temp
  # } /q

  eith <- 1 - foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                      .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
    {sapply(i, function(j) {CensCase__(data$exit[j], data$scr[[j]], 
                                          newpar[j, ], gkmat61$x, 
                                          gkmat61$w, glmat40$x, glmat40$w)})} / q
  # noscr <- 1 - foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                              .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
  #   {sapply(i, function(i) {CensCase__(data$exit[i], data$e.scr[[i]], 
  #                                     newpar[i, ], gkmat61$x, 
  #                                     gkmat61$w, glmat40$x, glmat40$w)})} / q
  int <- eith - sc
  stopCluster(cl)
  is_screen <- as.numeric(data$mode == 1)
  is_screen[is.na(is_screen)] <- 0
  is_interval <- as.numeric(data$mode == 0)
  is_interval[is.na(is_interval)] <- 0
  
  out <- cbind(data$studieid, data$entry, q, sapply(data$scr, max), #noscr, 
               eith, sc, int, 
               #(d10_int1 + d10_int2noscr), #/noscr, 
               #(d10_scr + d10_int1 + d10_int2scr), #/eith, 
               data$case, is_screen, is_interval)
  #colnames(out) <- c('studieid', 'age', 'surv', 'next_scr', 'risk_noscr', 'risk_overall', 'risk_screen', 'risk_interval', 'd10_noscr', 'd10_scr', 'case', 'is_screen', 'is_interval')
  cat("Prediction successful!\n")
  list(pred = as.data.frame(out), scr_size = matrix(screensize, length(ix), byrow = T), int_size = matrix(intsize1, length(ix), byrow = T) + matrix(intsize2, length(ix), byrow = T))
  # as.data.frame(out)
}

SizePred <- function(model, data, d_max = 10,  #x, c, s, scr, par,
                     gk1, gl, gk2 = NULL, followup = NULL,
                     d0 = 0.5, v0 = 0.06544985, n.cores = parallel::detectCores() - 1){
  
  library(foreach)
  library(doParallel)
  library(itertools)
  
  if(is.null(gk2)){
    gk2 <- gk1
  }
  if(is.null(followup)) {
    followup <- data$exit
  } else {
    followup <- data$entry + followup
  }
  
  v_max <- dv(d_max)
  v0 <- dv(d0)
  
  # cat("Predicting...\n")
  newpar <- BCmodelParBuilder(model, data)
  ix <- seq_along(data$entry)
  cl <- makePSOCKcluster(n.cores)
  registerDoParallel(cl)
  q <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
               .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
    sapply(i, function(i) {BC::f_CensCase(data$entry[i],data$e.scr[[i]], newpar[i, ],
                                      gk1$x, gk1$w, gl$x, gl$w)})
  d10_int1 <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                       .packages = "BC", .export = c("dv")) %dopar% {
                         sapply(i, function(j) {
                           ud.grid <- expand.grid(u=gk2$x * (data$scr[[j]][1]-data$entry[j])/2 + (data$scr[[j]][1]+data$entry[j])/2,
                                                  d=gk2$x * (d_max - d0)/2 + (d_max + d0)/2)
                           wts <- (gk2$w * (d_max-d0)/2) %*% t((data$scr[[j]][1]-data$entry[j])/2 * gk2$w)
                           temp <- mapply(function(u, d) {BC::f_SymptCase(u, dv(d), data$e.scr[[j]], newpar[j, ], gk1$x, gk1$w) * pi * d^2 / 2 },
                                          ud.grid$u, ud.grid$d)
                           sum(matrix(temp, nrow(wts), , T) * wts)
                         })}  / q
  d10_int2scr <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                          .packages = "BC", .export = c("dv")) %dopar% {
                            sapply(i, function(j) {
                              ud.grid <- expand.grid(u=gk2$x * (followup[j]-data$scr[[j]][1])/2 + (data$scr[[j]][1]+followup[j])/2,
                                                     d=gk2$x * (d_max - d0)/2 + (d_max + d0)/2)
                              wts <- (gk2$w * (d_max-d0)/2) %*% t((followup[j]-data$scr[[j]][1])/2 * gk2$w)
                              temp <- mapply(function(u, d) {BC::f_SymptCase(u, dv(d), data$scr[[j]], newpar[j, ], gk1$x, gk1$w) * pi * d^2 / 2 },
                                             ud.grid$u, ud.grid$d)
                              sum(matrix(temp, nrow(wts), , T) * wts)
                            })}  / q
  d10_int2noscr <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                            .packages = "BC", .export = c("gkmat61", "gkmat11", "dv")) %dopar% {
                              sapply(i, function(j) {
                                ud.grid <- expand.grid(u=gk2$x * (followup[j]-data$scr[[j]][1])/2 + (data$scr[[j]][1]+followup[j])/2,
                                                       d=gk2$x * (d_max - d0)/2 + (d_max + d0)/2)
                                wts <- (gk2$w * (d_max-d0)/2) %*% t((followup[j]-data$scr[[j]][1])/2 * gk2$w)
                                temp <- mapply(function(u, d) {BC::f_SymptCase(u, dv(d), data$e.scr[[j]], newpar[j, ], gk1$x, gk1$w) * pi * d^2 / 2 },
                                               ud.grid$u, ud.grid$d)
                                sum(matrix(temp, nrow(wts), , T) * wts)
                              })}  / q
  d10_scr <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                     .packages = "BC", .export = c("gkmat61", "gkmat11", "dv")) %dopar% {
                       dj.grid <- expand.grid(d=gk2$x * (d_max - d0)/2 + (d_max + d0)/2,
                                              j=i)
                       temp <- mapply(function(d, j) {pi * d^2 / 2 * BC::f_ScreenCase(data$scr[[j]][1], dv(d), data$scr[[j]], newpar[j, ], gk1$x, gk1$w)},
                                      dj.grid$d, dj.grid$j)
                       matrix(temp, , nrow(gk2), byrow = T) %*% (gk2$w * (d_max - d0)/2)
                       #temp
                     } /q
  # d10_int1 <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                      .packages = "BC", .export = c("d0", "dv")) %dopar% {
  #                        sapply(i, function(j) {
  #                          ud.grid <- expand.grid(u=gk2$x * (data$scr[[j]][1]-data$entry[j])/2 + (data$scr[[j]][1]+data$entry[j])/2,
  #                                                 d=gk2$x * (v_max - v0)/2 + (v_max + v0)/2)
  #                          wts <- (gk2$w * (v_max-v0)/2) %*% t((data$scr[[j]][1]-data$entry[j])/2 * gk2$w)
  #                          temp <- mapply(function(u, d) {BC::f_SymptCase(u, d, data$e.scr[[j]], newpar[j, ], gk1$x, gk1$w) },
  #                                         ud.grid$u, ud.grid$d)
  #                          sum(matrix(temp, nrow(wts), , T) * wts)
  #                        })}  / q
  # d10_int2scr <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                         .packages = "BC", .export = c("d0", "dv")) %dopar% {
  #                           sapply(i, function(j) {
  #                             ud.grid <- expand.grid(u=gk2$x * (data$exit[j]-data$scr[[j]][1])/2 + (data$scr[[j]][1]+data$exit[j])/2,
  #                                                    d=gk2$x * (v_max - v0)/2 + (v_max + v0)/2)
  #                             wts <- (gk2$w * (v_max-v0)/2) %*% t((data$exit[j]-data$scr[[j]][1])/2 * gk2$w)
  #                             temp <- mapply(function(u, d) {BC::f_SymptCase(u, d, data$scr[[j]], newpar[j, ], gk1$x, gk1$w) },
  #                                            ud.grid$u, ud.grid$d)
  #                             sum(matrix(temp, nrow(wts), , T) * wts)
  #                           })}  / q
  # d10_int2noscr <-  foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                           .packages = "BC", .export = c("gkmat61", "gkmat11", "d0", "dv")) %dopar% {
  #                             sapply(i, function(j) {
  #                               ud.grid <- expand.grid(u=gk2$x * (data$exit[j]-data$scr[[j]][1])/2 + (data$scr[[j]][1]+data$exit[j])/2,
  #                                                      d=gk2$x * (v_max - v0)/2 + (v_max + v0)/2)
  #                               wts <- (gk2$w * (v_max-v0)/2) %*% t((data$exit[j]-data$scr[[j]][1])/2 * gk2$w)
  #                               temp <- mapply(function(u, d) {BC::f_SymptCase(u, d, data$e.scr[[j]], newpar[j, ], gk1$x, gk1$w)},
  #                                              ud.grid$u, ud.grid$d)
  #                               sum(matrix(temp, nrow(wts), , T) * wts)
  #                             })}  / q
  # d10_scr <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                    .packages = "BC", .export = c("gkmat61", "gkmat11", "d0", "dv")) %dopar% {
  #                      dj.grid <- expand.grid(d=gk2$x * (v_max - v0)/2 + (v_max + v0)/2,
  #                                             j=i)
  #                      temp <- mapply(function(d, j) {BC::f_ScreenCase(data$scr[[j]][1], d, data$scr[[j]], newpar[j, ], gk1$x, gk1$w)},
  #                                     dj.grid$d, dj.grid$j)
  #                      matrix(temp, , nrow(gk2), byrow = T) %*% (gk2$w * (v_max - v0)/2)
  #                      #temp
  #                    } /q
  
  overall <- 1 - foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                     .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar% {
    sapply(i, function(i) {BC::f_CensCase(followup[i], data$scr[[i]], newpar[i, ],
                                      gkmat61$x, gkmat61$w, glmat40$x, glmat40$w) })} / q
  no_overall <- 1 - foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                         .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar% {
                           sapply(i, function(i) {BC::f_CensCase(followup[i], data$e.scr[[i]], newpar[i, ],
                                                                 gkmat61$x, gkmat61$w, glmat40$x, glmat40$w) })} / q
  # attend = foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                  .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
  #   sapply(i, function(i) {BC::f_SizeBounded(data$exit[i], v_max, data$entry[i], data$newscr[[i]], data$e.scr[[i]], newpar[i, ],
  #                                           gkmat61$x, gkmat61$w, glmat40$x, glmat40$w) / q[i]})
  # noattend = foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
  #                  .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
  #   sapply(i, function(i) {BC::f_SizeBounded(data$exit[i], v_max, data$entry[i], 255, data$e.scr[[i]], newpar[i, ],
  #                                            gkmat61$x, gkmat61$w, glmat40$x, glmat40$w) / q[i]})
  stopCluster(cl)
  # cat("Prediction successful!\n")
  data.frame(attend = d10_int1 + d10_scr + d10_int2scr, no_attend = d10_int1 + d10_int2noscr, 
             overall = overall, no_overall = no_overall, entry = q)
  # cbind(d10_int1 + d10_int2noscr, d10_int1 + d10_scr + d10_int2scr)
}

BCmodelLLtest <- function(model1, model0) {
  if (!(model1$fitted && model0$fitted)){
    print("At least one model has not been fitted!")
    out <- NULL
  } else {
    t <- 2 * abs(model1$LogL - model0$LogL)
    k <- abs(length(model1$par) - length(model0$par))
    out <- data.frame("statistic" = round(t, 2), "df" = k, "p" = round(1 - pchisq(t, k), 5))
  }
  out
}

BCsize <- function(model, data, n.cores = parallel::detectCores() - 1) {
  library(foreach)
  library(doParallel)
  library(itertools)
  
  RiskScreen__ <- ifelse(model$OnsetDependentGrowth, BC::f_RiskScreen_odg, BC::f_RiskScreen)
  CensCase__ <- ifelse(model$OnsetDependentGrowth, BC::f_CensCase_odg, BC::f_CensCase)
  ScreenCase__ <- ifelse(model$OnsetDependentGrowth, BC::f_ScreenCase_odg, BC::f_ScreenCase)
  SymptCase__ <- ifelse(model$OnsetDependentGrowth, BC::f_SymptCase_odg, BC::f_SymptCase)
  
  cat("Predicting...\n")
  newpar <- BCmodelParBuilder(model, data)
  ix <- seq_along(data$entry)
  cl <- makePSOCKcluster(n.cores)
  registerDoParallel(cl)
  q <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
               .packages = "BC", .export = c("gkmat61", "glmat40")) %dopar%
    sapply(i, function(i) {CensCase__(data$entry[i],data$e.scr[[i]], newpar[i, ],
                                      gkmat61$x, gkmat61$w, glmat40$x, glmat40$w)})
  screensize <- foreach(i=isplitVector(ix, chunks = n.cores), .combine = c,
                        .packages = "BC", .export = c("gkmat61", "glmat40", "d0", "dv")) %dopar% {
                          jv.grid <- expand.grid(v=glmat40$x, j=i)
                          mapply(function(j, d) {pi * (d+d0)^2 / 2 * ScreenCase__(data$entry[j], dv(d + d0), c(data$entry[j], data$e.scr[[j]]), newpar[j, ], gkmat61$x, gkmat61$w) / q[j]},
                                 jv.grid$j, jv.grid$v)}
  #matrix(temp, length(i), , F)
  

  stopCluster(cl)
  
  cat("Prediction successful!\n")
  matrix(screensize, length(ix), byrow = T)
  #as.data.frame(out)
}

BCprofile <- function(pdcat, famhist = 0, relhist = 0, bmi_m = 0, 
                      schedule = seq(40, 74, by = 2), del = 1/52) {
  
  schedule_ex <- c(0, schedule, 80)
  age.cat <- cut(p3data$entry, c(0, 2 * (20:38) - 1, 100))
  
  pd.q <- NA
  pd.q <- t(sapply(levels(age.cat), function(each) {
    ix <- which(age.cat == each)
    pd.lvl <- p3data$pd[ix] #+ rnorm(length(ix), 0, 1e-5)
    quantile(pd.lvl, 0.01 * c(25, 50, 75, 100) - 0.125, na.rm = T)
  }))
  
  df_full <- NULL
  for(k in 2:length(schedule_ex)){
    age <- seq(schedule_ex[k-1] + del, schedule_ex[k], by = del)
    scr <- lapply(1:length(age), function(i) schedule[which(schedule < schedule_ex[k])])
    df <- data.frame(age = age)
    df$scr <- scr
    df_full <- rbind(df_full, df) 
  }
  
  df_full$agecat <- cut(df_full$age, c(0, 2 * (20:38) - 1, 100))
  df_full$agecat[1] <- df_full$agecat[2] 
  df_full$pd1 <- pd.q[df_full$agecat, 1]
  df_full$pd2 <- pd.q[df_full$agecat, 2]
  df_full$pd3 <- pd.q[df_full$agecat, 3]
  df_full$pd4 <- pd.q[df_full$agecat, 4]
  df_full$pdcat <- pdcat
  df_full$pdcat2 <- df_full$pdcat == 2
  df_full$pdcat3 <- df_full$pdcat == 3
  df_full$pdcat4 <- df_full$pdcat == 4
  df_full$pd <- pd.q[df_full$agecat, df_full$pdcat[1]]
  df_full$famhist  <- famhist
  df_full$relhist  <- relhist
  df_full$bmi_m  <- bmi_m
  
  newpar <- BCmodelParBuilder(mymodel7, df_full)
  test <- sapply(1:nrow(df_full), function(i) {
    BC::f_CensCase(df_full$age[i], df_full$scr[[i]],
                   newpar[i, ], gkmat61$x, gkmat61$w, glmat40$x, glmat40$w)
  })
  test2 <- sapply(1:nrow(df_full), function(i) {
    BC::f_CensCase(df_full$age[i], df_full$scr[[1]],
                   newpar[i, ], gkmat61$x, gkmat61$w, glmat40$x, glmat40$w)
  })
  cbind(age = df_full$age, screened = test, unscreened = test2)
}

