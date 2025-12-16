#### libraries
wd=getwd()
source(paste0(wd, "/Methods/helper/loader.R"))
source(paste0(wd, "/Methods/helper/wrapper_point.R"))
source(paste0(wd, "/Methods/helper/wrapper_coverage.R"))
source(paste0(wd, "/Methods/helper/tabler.R"))
source(paste0(wd, "/Methods/helper/binder.R"))
source(paste0(wd, "/Methods/helper/melter.R"))
source(paste0(wd, "/Methods/helper/postestimation.R"))

library(kableExtra)


#### Utility vectors. to be used inside the helper functions to estimate
# test statistics for different time lenghts

tt0=c(rep(10,3), rep(20, 3), rep(40,3))
pt1=rep(5,9)
pt2=rep(10,9)
pt3=c(rep(5,3), rep(10, 3), rep(20,3))

methods=c("SC", "SR", "BVR", "BSC", "SMAC")
scenarios <- c("t0=10, IID", "t0=10, 40%","t0=10, 70%",
               "t0=20, IID", "t0=20, 40%","t0=20, 70%",
               "t0=40, IID", "t0=40, 40%","t0=40, 70%")
wd<-"C:/Users/giuli/Desktop/cluster ottobre"
## -- POINT ESTIMATES RESULTS -- ##
point_001=postestimation(wd,  "0.01",   mode = 1,  scenarios)
point_02=postestimation(wd,  "0.2",   mode = 1,  scenarios)
point_04=postestimation(wd,  "0.4",   mode = 1,  scenarios)
point_06=postestimation(wd,  "0.6",   mode = 1,  scenarios)

## -- COVERAGE -- ##

cov_04=postestimation(wd,  "0.4",   mode = 2,  scenarios)
cov_001=postestimation(wd,  "0.01",   mode = 2,  scenarios)
cov_02=postestimation(wd,  "0.2",   mode = 2,  scenarios)
cov_06=postestimation(wd,  "0.6",   mode = 2,  scenarios)

## TABLER FUNCTION - it generates the table to be used in kable
################################################################################
## we need to generate three different tables according to the length
## of post-treatment
#################################

## bias tables

tab_bias001<- cbind(tabler(point_001, tt0, pt1, stat="bias"), 
                   tabler(point_001, tt0, pt2, stat="bias"),
                   tabler(point_001, tt0, pt3, stat="bias"))  

kable(round(tab_bias001,3), format="latex",booktabs=T)

tab_bias04<- cbind(tabler(point_04, tt0, pt1, stat="bias"), 
                   tabler(point_04, tt0, pt2, stat="bias"),
                   tabler(point_04, tt0, pt3, stat="bias"))  

kable(round(tab_bias04,3), format="latex",booktabs=T)

tab_bias06<- cbind(tabler(point_06, tt0, pt1, stat="bias"), 
                   tabler(point_06, tt0, pt2, stat="bias"),
                   tabler(point_06, tt0, pt3, stat="bias"))  

kable(round(tab_bias06,3), format="latex",booktabs=T)

## MSE tables

tab_mse001<- cbind(tabler(point_001, tt0, pt1, stat="MSE"), 
                  tabler(point_001, tt0, pt2, stat="MSE"),
                  tabler(point_001, tt0, pt3, stat="MSE"))  

kable(round(tab_mse001,3), format="latex",booktabs=T)

tab_mse04<- cbind(tabler(point_04, tt0, pt1, stat="MSE"), 
                   tabler(point_04, tt0, pt2, stat="MSE"),
                   tabler(point_04, tt0, pt3, stat="MSE"))  

kable(round(tab_mse04,3), format="latex",booktabs=T)

tab_mse06<- cbind(tabler(point_06, tt0, pt1, stat="MSE"), 
                  tabler(point_06, tt0, pt2, stat="MSE"),
                  tabler(point_06, tt0, pt3, stat="MSE"))  

kable(round(tab_mse06,3), format="latex",booktabs=T)



## Coverage tables
tab_cov001<- cbind(tabler(cov_001, tt0, pt1, stat="coverage"),
                  tabler(cov_001, tt0, pt2, stat="coverage"),
                  tabler(cov_001, tt0, pt3, stat="coverage"))  

kable(round(tab_cov001,3), format="latex",booktabs=T)

tab_cov02<- cbind(tabler(cov_02, tt0, pt1, stat="coverage"),
                  tabler(cov_02, tt0, pt2, stat="coverage"),
                  tabler(cov_02, tt0, pt3, stat="coverage"))  

kable(round(tab_cov02,3), format="latex",booktabs=T)

tab_cov04<- cbind(tabler(cov_04, tt0, pt1, stat="coverage"),
                  tabler(cov_04, tt0, pt2, stat="coverage"),
                  tabler(cov_04, tt0, pt3, stat="coverage"))  

kable(round(tab_cov06,3), format="latex",booktabs=T)


tab_cov06<- cbind(tabler(cov_06, tt0, pt1, stat="coverage"),
                  tabler(cov_06, tt0, pt2, stat="coverage"),
                   tabler(cov_06, tt0, pt3, stat="coverage"))  

kable(round(tab_cov06,3), format="latex",booktabs=T)

### end of tables
