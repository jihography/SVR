######################################################
########################################################
## master script for computations on hipergator  #######
########################################################
########################################################

args <- commandArgs(TRUE)
args
sp_range <- as.numeric(args[1])
tt_periods <- as.numeric(args[2])
errors_sp <- as.numeric(args[3])
index <- as.numeric(args[4])
set.seed(index)

setwd("/Users/kjh/Documents/innocity/SVR")

# Set to TRUE to bypass simulations and load an external dataset that already
# contains a matrix named `sim` (time x units). You can optionally include
# `bands`, `t0`, and `num_controls` in the .RData file to override the defaults
# below. The path is relative to the working directory set above.
use_external_data <- FALSE
external_data_path <- "SVR_input_data.RData"

# --------------- Sourcing functions --------------- #

# Sourcing in code from other files.
source("Functions/sim_model_function.R")
source("Functions/estimation_function.R")
source("Functions/calculation_function.R")
source("Functions/point_estimate_function.R")
# source("Functions/postestimation_function.R")  # Not found.
source("Functions/ci_function.R")
source("Functions/coverage_function.R")
source("Methods/helper/preset_function.R")
source("Methods/helper/ci_shen.R")
source("Methods/helper/ci_bayes.R")
source("Methods/helper/cover_helper.R")
source("Methods/helper/binder.R")
source("Methods/helper/wrapper.R")
source("Methods/helper/helper_vertical_regression.R")

# Function for performing synthetic controls.
source("Methods/helper/SCM_function.R")
source('Methods/sepSC_method.R')
source('Functions/helper/sepSC_calc_function.R')

# Function for performing separate vertical regression with ridge.
source('Methods/sepSR_method.R')
source('Functions/helper/sepSR_calc_function.R')

# Function for performing separate vertical regression with SVD
source('Methods/sepOLS_method.R')

# Function for performing Bayesian separate vertical regression.
source('Methods/sepBVR_method.R')

# Function for performing Bayesian synthetic control.
source('Methods/sepBSC_method.R')

# Function for performing SMAC. 
source('Methods/SMAC_method.R')

# ------------ End of Sourcing functions ------------- #


library(LowRankQP)
library(glmnet)
library(rstan)
library(spatstat)
library(Matrix)
library(fungible)
rstan_options(auto_write = FALSE)
#out_path <- 'Output/1_sims/Results/'


# ----------- PART A: Setting the simulation parameters ----------- #

## dim parameter
num_controls <- 7
t0 <- tt_periods
time_periods <- t0 + 20
time_periods_controls <- 80  # -GP- Do not change this with t0.
bands <- 2

## sampling pars
iter <- 6000
warm <- 2000
chains <- 3
## sim pars

# ----- Spatial and temporal correlation parameters.

sp_var <- .4
tt_var <- 0.3 ^ 2  # -GP- I reduced the tt_var a little.
ti_var <- 0.7 ^ 2  # -GP- I squared this to make it a variance and comparable to tt_var
bi_var <- 0.5 ^ 2  # -GP- I squared this to make it a variance and comparable to sp_var
tt_range <- .05
sp_nugget <- 0.001
tt_nugget <- 0.15 ^ 2  # -GP- Adding some temporal nugget to the controls.
rho_error <- .2	  


# ----- Outcome model errors.

if (errors_sp == 1) {
  e_weight <- 0  # Proportion of error's variance that is spatial
  share_error <- 0.4  # Noise-signal ratio in terms of variances (error sd as % of signal)
} else if (errors_sp == 2) {
  e_weight <- .5
  share_error <- 0.4
} else if (errors_sp == 3) {
  e_weight <- .5
  share_error <- 0.7
}
print(errors_sp)


# ------- The methods that will be used.
method <- c("SC","SR", "OLS", "BVR", "BSC", "SMAC")
method <- c("SC","SR", "OLS")

# --------- The treated units.
treated_radius <- seq(0, 1, length.out = bands)
print(treated_radius)


# ---------------- PART B: Generating or loading data ---------------- #

if (use_external_data) {
  external_env <- new.env()
  loaded_vars <- load(external_data_path, envir = external_env)

  if (!"sim" %in% loaded_vars) {
    stop("The external data file must contain an object named 'sim'.")
  }

  sim <- external_env$sim

  if ("bands" %in% loaded_vars) {
    bands <- external_env$bands
  }

  if ("t0" %in% loaded_vars) {
    t0 <- external_env$t0
    tt_periods <- t0
  }

  if ("num_controls" %in% loaded_vars) {
    num_controls <- external_env$num_controls
  } else {
    num_controls <- ncol(sim) - bands
  }

  time_periods <- nrow(sim)
  beta_true <- NULL

  message("Loaded external sim matrix with ", nrow(sim), " rows and ",
          ncol(sim), " columns.")
} else {
  seed_b <- seed_t <- seed_e <- index

  sim <- sim_model(seed_b = seed_b, seed_t = seed_t, seed_e = seed_e,
                   time_periods = time_periods,
                   time_periods_controls = time_periods_controls,
                   bands = bands, num_controls = num_controls,
                   sp_var = sp_var, sp_range = sp_range, bi_var = bi_var,
                   tt_var = tt_var,
                   tt_range = tt_range, ti_var = ti_var, sp_nugget = sp_nugget,
                   tt_nugget = tt_nugget,
                   e_weight = e_weight, share_error = share_error)

  beta_true <- sim$beta
  sim <- sim$sim  # The potential outcomes under control.

  set.seed(index)
}



# ---------------- PART C: Estimating the models ---------------- #

iter <- 4
warm <- 2
chains <- 1

est <- estimation(sim = sim, t0 = t0, bands = bands, iter = iter, warm = warm,
                  norm = TRUE, method = method, chains = chains)


# ---------------- PART D: Getting predictions ---------------- #

cal <- calculation(sim = sim, est = est, bands = bands, norm = TRUE)
point <- point_estimate(sim, cal) ## it calculates bias and MSE
c_interv <- ci(sim = sim, est = est, cal = cal, t0 = t0, norm = TRUE)
cover <- coverage(sim, interv = c_interv)



# ---------------- PART E: Saving results ---------------- #

res <- list(sim = sim, est = est, cal = cal, beta_true = beta_true,
            point = point, ci = c_interv, coverage = cover)

out_path <- paste0('Output/apr_sims/Results/ss', sp_range, '/tt', tt_periods,
                   '/ee', errors_sp)
out_path
out_file=paste0(out_path, "/res_",index, ".RData")
out_file
save(res, file=out_file)
print("fine primo esperimento")

################################################################################

