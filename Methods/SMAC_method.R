#' Performing separate SMAC. 
#' 
#' @param ym.pre Matrix T0 x bands for the pre-intervention period and the
#' outcomes units.
#' @param x.pre Matrix T0 x (# controls) for the pre-intervention period and
#' the control units.
#' @param x Matrix T x (# controls) for the pre-intervention and 
#' post-intervention period and the control units.
#' @param treated_radius vector of (# treated units) distances across treated
#' @param chains How many MCMC chains to run. Defaults to 3.
#' 
#' Depends on the MGP.stan file, and rstan library
#' 
SMAC <- function(ym.pre, x.pre, x, treated_radius, chains = 3) {
  
  # arguments
  bands <- ncol(ym.pre)
  time_periods <- nrow(x)
  t0 <- nrow(x.pre)
  
  # Without an intercept:
  ss_data = list(
    x = x.pre,
    y = ym.pre,
    xnn = x,
    h = treated_radius,
    tt = time_periods,
    t0 = t0,
    n_t = bands, 
    n_c = num_controls
  )
  
  # With an intercept it would be:
  #   x = cbind(1, x.pre),
  #   xnn = cbind(1, x),
  #   n_c=num_controls+1
  
  fit <- rstan::stan(
    file = "Methods/SMAC.stan",  # SMAC.
    data = ss_data,
    cores = min(chains, 3),
    iter = iter,
    chains = chains,
    verbose = F,
    warmup = warm,
    control = list(
      max_treedepth = 15,
      stepsize = 0.02,
      adapt_delta = 0.99
    )
  )
  
  parameters <- rstan::extract(fit)
  
  return(parameters)
  
}