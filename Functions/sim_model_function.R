################################################################################
################## MODEL BASED SIMULATIONS #####################################
################################################################################



### inputs 
#' @sp_var=1  ## marginal spatial variance
#' @sp_range=1 ## marginal spatial lengthscale - rho
#' @tt_var=1 ## marginal temporal variance 
#' @tt_range=.01 ##marginal temporal lengthscale - rho
#' @sp_nugget=.1 ## spatial nugget
#' @tt_nugget=.1 ## temporal nugget
### time periods - number of treated units - number of controls 
### no need to specify the number of t0 because all the process should be GP


sim_model<-function(seed_b, seed_t, seed_e,
                    time_periods, time_periods_controls, bands, num_controls,
                    treated_radius, rho_error,
                    sp_var, sp_range, bi_var,  tt_var, tt_range, ti_var, sp_nugget, tt_nugget,
                    e_weight, share_error){
  
  
  # ----------- PART A: Creating the control units across time. ------------ #
  
  ## construct a matrix of control units from a GP - BUT COULD BE ALSO SOMETHING DIFFERENT  
  y = matrix(NA, nrow = time_periods_controls, ncol = num_controls)
  dimnames(y) <- list(time = 1 : time_periods_controls, control = 1 : num_controls)
  set.seed(seed_t)
  
  # Baseline values for the control time series. Each unit can vary around a
  # different overall value.
  ti <- rnorm(num_controls, 0, sqrt(ti_var))
  
  # Creating the temporal covariance matrix of the controls across time
  # -GP- Lines 36-41 were previously inside the loop. Putting them outside will speed things up.
  curr_alpha <- tt_var
  curr_rho <- tt_range
  xpoints <- seq(0, 1, length.out = time_periods_controls)      ## distance to be used - temporal GP 
  xpoints<-xpoints/max(xpoints)
  D_xloc <- as.matrix(dist((xpoints)))                   ## temporal GP 
  Sigma <- (curr_alpha * (exp(- D_xloc ^ 2 / (2*(curr_rho ^ 2)))) +
              tt_nugget * diag(1, (time_periods_controls))) ## covariance matrix of the GP
  
  for (bb in 1 : num_controls) {
    y[, bb] <- ti[bb] + t(chol(Sigma)) %*% as.vector(rnorm(time_periods_controls, 0, 1)) ## consruction of the control time series
  }
  
  # Truncating to only the time periods we want to keep:
  y <- y[1 : time_periods, ]
  
  
  # ---- PART B: Creating the coefficients of the controls across space ----- #
  
  #### construct a matrix of GP coefficients.
  beta <- matrix(nrow = bands, ncol = num_controls)
  dimnames(beta) <- list(treated = 1 : bands, control = 1 : num_controls)
  set.seed(seed_b)
  
  # bi below represents an overall importance for control i. This overall
  # importance can change across space according to a GP.
  bi <- rnorm(num_controls, 0, sqrt(bi_var))
  
  # Creating the covariance matrix of the coefficients across space.
  # -GP- This does not need to be in the loop.
  curr_alpha <- sp_var  # fixed variance
  curr_rho <- sp_range   # fixed lengthscale
  D_xloc <- as.matrix(dist(sort(treated_radius))) # distance to be used - spatial GP 
  # covariance matrix of the beta GP 
  Sigma <- (curr_alpha * (exp(- D_xloc ^ 2 / (2*(curr_rho ^ 2)))) + 
              sp_nugget * diag(1, length(treated_radius)))
  
  for (bb in 1 : num_controls) {
    beta[, bb] <- bi[bb] + t(chol(Sigma)) %*% as.vector(rnorm(bands, 0,1)) # construction of beta GP 
  }
  
  
  # ---- PART C: Generating potential outcomes for the treated units ----- #
  
  # (A) The potential outcome under control:
  
  y_tr = y %*% t(beta)
  dimnames(y_tr) <- list(time = 1 : time_periods, treated = 1 : bands)
  
  true_sim <- cbind(y_tr, y)  # Dimension T x (number of units, treated first)
  colnames(true_sim) <- c(paste('treated', 1 : bands), paste('control', 1 : num_controls))
  names(dimnames(true_sim)) <- c('time', 'unit')
  
  # (B) Adding the error terms for the treated units.
  
  # Getting the variance of the treated units to get the "signal"
  var_true <- mean(apply(true_sim[, 1 : bands], 2, var))
  # Proportion of signal for the standard deviation error term:
  error <- sqrt(share_error * var_true)
  
  ### errors
  ## from 27/11 it does not count anymore the boolean for spatial but the weight we are assigning to the spatial part 
  
  # The spatial error for the treated units (variance 1).
  set.seed(seed_e) # -GG- : set the seed
  errors_1 <- matrix(NA, nrow = time_periods, ncol = bands)
  for (i in 1 : time_periods) {
    Sigma <- (exp(- D_xloc ^ 2 / (2*(rho_error ^ 2)))) + sp_nugget * diag(1, bands)
    errors_1[i,] <- t(chol(Sigma)) %*% as.vector(rnorm(bands, 0, error)) # construction of beta GP
  }
  errors_0 <- matrix(0, nrow = time_periods, ncol = num_controls)
  # errors includes the spatial error term for the treated units and 0s for the
  # control units.
  errors <- cbind(errors_1, errors_0)
  colnames(errors) <- c(paste('treated', 1 : bands), paste('control', 1 : num_controls))
  
  # The non-spatial error for the treated units (variance 1).
  wn <- rnorm(bands * time_periods, 0, error)
  wn <- matrix(wn, nrow = time_periods)
  # wn includes the non-spatial error term for the treated units and 0s for the
  # control units.
  wn <- cbind(wn, errors_0)
  colnames(wn) <- colnames(errors)
  
  # - GP! - Don't change the line below.
  ee <- sqrt(e_weight) * errors + sqrt(1 - e_weight) * wn
  #corrplot::corrplot(cor(errors))
  
  sim <- true_sim + ee
  #corrplot::corrplot(cor(sim))  # check visual results
  
  out <- list()
  out[[1]] <- sim  # The potential outcomes (under control) for treated and control units.
  out[[2]] <- beta  # The true coefficients of the controls across radii.
  out[[3]] <- y  # The observed control time series. (Also part of sim.)
  out[[4]] <- true_sim  # The underlying signal of the time series, w/o error added to treated units (unobservable).
  names(out) <- c('sim', 'beta', 'y', 'true_sim')
  return(out)
}


