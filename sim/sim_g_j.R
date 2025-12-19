

library(Rcpp)
library(eha)
library(flexsurv)
library(parallel)
sourceCpp("./utility/commonf.cpp")
source("./utility/dt_gen_fct.r")
source("./sim/utility/dt_csd_gen.r")
source("./sim/utility/commonf.r")
source("./utility/draft_est_rec20_g_j.r")

n <- 10000
#lognormal parameters for age of first sex
meanlog = 2.85
sdlog = 0.15
d <- 0 #duration of natural immunity
#exponential parameters for appearance
cuts=22
lam01=c(0.2, 0.05)
#q01_22 <- .2; q01_22_29 <- .05 ;q01_30p <- .05

#weibull parameters for time from appearance to clearance
q10_scale <- 0.75
q10_shape <- 0.74
#weibull parameters for progression to precancer
q12_scale <- 4.75#15/2
q12_shape <- 2.2#1.2
#weibull parameters for regression of precancers
q20_scale <- 10#2.5#2.5
q20_shape <- 1#0.4

sim.fun <- function(n.sim, output1,seed){
  
  set.seed(seed)
  
  for(k in 1:n.sim){
    
    #set.seed(k)
    # data generation
    dat=dt_csd_gen.f(n, meanlog, sdlog, cuts, lam01, q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape)
    dtw<-aggregate(w~age+x,dat,sum)
    
    #dtw
    # estimation
    ini=c(rep(log(0.01),3))
    est=est_lam01_pwc01_wei_rec20_g_j.f(ini, dtw, cuts=c(22), lam12=c(q12_scale, q12_shape), lam10=c(q10_scale, q10_shape), meanlog, sdlog, Q=0, sexactive=1)
    
    saveRDS(c(est$est$est), file=output1)
    
    print(k)
  }
  
}
