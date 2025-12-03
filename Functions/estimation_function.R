###############################################################################
######################### COEFFICIENT ESTIMATION ##############################
###############################################################################

## @ sim <- data or simulation matrix
## @ est <- results from estimation function
## @ norm <- normalization of outcomes - pre treatment


estimation <- function(sim, t0, bands, treated_radius, iter, warm,
                       norm = TRUE, method, chains = 3) {
  
  # Getting some quantities that are used later in the code.
  num_controls <- dim(sim)[2] - bands
  
  if (norm) {  # Standardizing using the mean and SD in the pre-intervention period.
    sim <- preset(sim, t0)$std
  }
  
  # Preparing different matrix config for estimation:
  # Treated and control units, pre/post treatment, training
  
  # For the vectorized outcomes, all data for unit 1 are first, then data from
  # unit 2, etc. So as.vector loops over times and then units.
  
  # Getting the treated unit data:
  yv = as.vector(sim[, 1 : bands])  # As a vector over all time periods
  yv.pre = as.vector(sim[1 : t0, 1 : bands])  # As a vector in pre-intervention
  ym = sim[, 1 : bands]  # As a matrix over all time periods.
  ym.pre = sim[1 : t0, 1 : bands] # As a matrix in pre-intervention.
  
  # Getting the control unit data.
  x = sim[, (bands + 1) : (num_controls + bands)]  # Matrix, all time periods
  x.pre = sim[1 : t0, ((bands) + 1) : (num_controls + bands)]  # Matrix, pre-intervention
  
  # Getting the training data, will be used in ridge only.
  train = round(t0 - t0 / 5)
  y.train = as.vector(sim[1 : train, 1 : bands])
  x.train = sim[1:train, ((bands) + 1) : (num_controls + bands)]
  
  
  ### storing results
  
  out = list()
  
  ### 1.0 SEPARATED SCM
  
  if ("SC" %in% method){
    out$SC = sepSC(ym.pre = ym.pre, x.pre = x.pre)
    # The output is a matrix of dimension (num_controls x bands) with
    # estimated coefficients.
    print("SC estimates done")
  }
  
  ### 2- Separate ridge estimation
  
  if ("SR" %in% method){
    out$SR = sepSR(ym.pre = ym.pre, x.pre = x.pre)
    # The output is a matrix of dimension ((num_controls + 1) x bands) with
    # estimated coefficients. The first row is the intercept.
    print("SR estimates done")
  }
  
  if ("OLS" %in% method){ # Separate OLS for each treated unit.
    out$OLS = sepOLS(ym.pre = ym.pre, x.pre = x.pre)
    # The output is a matrix of dimension ((num_controls + 1) x bands) with
    # estimated coefficients. The first row is the intercept.
    print("OLS estimates done")
  }
  
  
  if ("BVR" %in% method) {
    out$BVR <- sepBVR(ym.pre = ym.pre, x.pre = x.pre, x = x, iter = iter,
                      warm = warm, chains = chains)
    # The output is a list, where each element of the list corresponds to a
    # treated unit. Then each element is a list of itself from what was
    # extracted from the Stan fit, including a vector of residual variances
    # (sigma_sq) and standard deviations (sigma), a matrix of coefficients
    # where the rows are iterations and the columns correspond to intercept
    # and controls (gp1), and the predicted values for the treated unit where
    # rows are iterations and columns are time periods (y_new).
    print("BVR estimates done")
  }
  
  
  if ("BSC" %in% method){
    out$BSC <- sepBSC(ym.pre = ym.pre, x.pre = x.pre, x = x, iter = iter,
                      warm = warm, chains = chains)
    # The exact same output as sepBVR, except the coefficient matrix does not
    # include intercepts.
    print("BSC estimates done")
  }
  
  
  if("SMAC" %in% method) {
    out$SMAC <- SMAC(ym.pre = ym.pre, x.pre = x.pre, x = x,
                     treated_radius = treated_radius, num_controls = num_controls,
                     iter = iter, warm = warm, chains = chains)
    # The output of the SMAC stan fit including coefficients, spatial
    # parameters for the errors and the coefficients, weight of spatial VS
    # iid error, covariance matrix, and predictions for the treated units.
    # The predictions are in ynn.
    print("SMAC estimates done")
  }
  
  
  library(purrr)
  out = compact(out)
  
  return(out)
  
}

