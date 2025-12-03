#' Performing separate Bayesian vertical regression for each band.
#' 
#' @param ym.pre Matrix T0 x bands for the pre-intervention period and the
#' outcomes units.
#' @param x.pre Matrix T0 x (# controls) for the pre-intervention period and
#' the control units.
#' @param x Matrix T x (# controls) for the pre-intervention and 
#' post-intervention period and the control units.
#' 
#' Depends on the BVR.stan file, and rstan library
#' 
sepBVR <- function(ym.pre, x.pre, x, iter, warm, chains = 3) {
  
  # arguments
  bands <- ncol(ym.pre)
  num_controls <- ncol(x.pre)
  parameters = list()

  bvr_model <- rstan::stan_model(file = "Methods/BVR.stan")
  has_parameters <- length(bvr_model@model_pars) > 0
  
  for (i in 1:bands){
    
    ss_data = list(
      x = cbind(1, x.pre),
      y = ym.pre[,i],
      N = nrow(ym.pre),
      C = num_controls + 1,
      X_new=cbind(1,x),
      N_new=nrow(x)
    )
    
    
    fit <- rstan::sampling(
      object = bvr_model,
      data = ss_data,
      cores = min(chains, 3),
      iter = iter,
      chains = chains,
      verbose = FALSE,
      warmup = warm,
      algorithm = if (has_parameters) "NUTS" else "Fixed_param",
      control = list(
        max_treedepth = 12,
        stepsize = 0.05,
        adapt_delta = 0.85
      )
    )
    parameters[[i]] <- rstan::extract(fit)
  }
  names(parameters) <- paste('treated', 1 : bands)
  return(parameters)
}

