###############################################################################
################ CALCULATION AFTER COEFFICIENT ESTIMATION #####################
###############################################################################
## @ sim <- data or simulation matrix
## @ est <- results from estimation function
## @ norm <- normalization of outcomes - pre treatment


calculation <- function(sim, est, t0, bands, norm) {
  
  # If the data are not standardized, we do not need to rescale.
  means <- rep(0, bands)
  sds <- rep(1, bands)
  
  # Standardizing the data using the preset function.
  if (norm) {
    sim_std <- preset(sim, t0)
    sim <- sim_std$std
    means <- sim_std$means[1 : bands]
    sds <- sim_std$sds[1 : bands]
  }
  
  # Preparing the control unit data.
  x = sim[, (bands + 1) : ncol(sim)]  # Matrix, all time periods
  
  
  out = list()
  
  ### 1.0 SEPARATED SCM
  if ("SC" %in% names(est)){
    out$SC = sepSC_calc(estSC = est[["SC"]], x = x, means = means, sds = sds) 
  }
  
  
  ### 2.0 SEPARATED VERTICAL REGRESSION WITH RIDGE PEN
  if ("SR" %in% names(est)){
    out$SR <- sepSR_calc(estSR = est[["SR"]], x = x, means = means, sds = sds)
  }
  
  ### 3.0 MULTIVARIATE OLS
  if ("OLS" %in% names(est)){
    out$OLS = x %*% as.matrix(est[["OLS"]])
    for (ii in 1 : bands) {
      out$OLS[, ii] <- out$OLS[, ii] * sds[ii] + means[ii]
    }
  }
  
  ### 4.0 BAYESIAN VERTICAL REGRESSION
  if ("BVR" %in% names(est)) {
    mat <- matrix(nrow = nrow(x), ncol = bands)  # For the posterior mean
    mat_med <- matrix(nrow = nrow(x), ncol = bands)  # For the posterior median
    for (i in 1 : bands) {
      samples <- est[["BVR"]][[i]]$y_new * sds[i] + means[i]
      mat[, i] <- apply(samples, 2, mean)
      mat_med[, i] <- apply(samples, 2, median)
    }
    out$BVR <- mat
    out$BVRmedian <- mat_med  # What we were keeping track before as BVR.
  }
  
  
  ### 5.0 BAYESIAN SYNTHETIC CONTROL
  if ("BSC" %in% names(est)) {
    mat <- matrix(nrow = nrow(x), ncol = bands)  # For the posterior mean
    mat_med <- matrix(nrow = nrow(x), ncol = bands)  # For the posterior median
    for (i in 1 : bands) {
      samples <- est[["BSC"]][[i]]$y_new * sds[i] + means[i]
      mat[, i] <- apply(samples, 2, mean)
      mat_med[, i] <- apply(samples, 2, median)
    }
    out$BSC <- mat
    out$BSCmedian <- mat_med  # What we were keeping track before as BSC.
  }
  
  
  ### 6.0 SMAC
  if ("SMAC" %in% names(est)){
    mat <- matrix(nrow = nrow(x), ncol = bands)  # For the posterior mean
    mat_med <- matrix(nrow = nrow(x), ncol = bands)  # For the posterior median
    for (i in 1 : bands) {
      samples <- est[["SMAC"]]$ynn[, , i] * sds[i] + means[i]
      mat[, i] <- apply(samples, 2, mean)
      mat_med[, i] <- apply(samples, 2, median)
    }
    out$SMAC <- mat
    out$SMACmedian <- mat_med
  }
  
  return(out)
}
