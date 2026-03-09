
library(MASS)
library(Rcpp)
library(eha)
library(flexsurv)
library(parallel)
sourceCpp("/home/maof3/HPV_appearance/sim1/utility/commonf.cpp")
source("/home/maof3/HPV_appearance/sim1/utility/dt_gen_fct.R")
source("/home/maof3/HPV_appearance/sim1/utility/commonf.R")
source("/home/maof3/HPV_appearance/sim1/utility/est_rec20.R")

n <- 5000
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
q12_scale <- 6.17#4.75#15/2
q12_shape <- 2.43#2.2#1.2
#weibull parameters for regression of precancers
q20_scale <- 12.5#2.5#2.5
q20_shape <- 1#0.4

sim.fun <- function(n.sim, output1,seed){
  
  set.seed(seed)
  
  for(k in 1:n.sim){
    
    #set.seed(k)
    # data generation
    dat=dt_csd_gen.f(n, meanlog, sdlog, cuts, lam01, q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape)
    dtw<-aggregate(w~age+hpv,dat,sum)

    # estimation
    ini=c(rep(log(0.01),2))
    #est=est01_20_g_S.f(ini=ini, dat, cuts=c(22), c(q10_scale, q10_shape), c(q12_scale, q12_shape), num.cores=2)
    est=est01.f(ini, dtw, cuts=c(22), lam12=c(q12_scale, q12_shape), lam10=c(q10_scale, q10_shape), lam20=c(q20_scale, q20_shape), 
		meanlog, sdlog, Q=0, sexactive=1)
    
    swvar=swvar_01.f(gradstep=1e-05, est$est, dtw, cuts=c(22), lam12=c(q12_scale, q12_shape), lam10=c(q10_scale, q10_shape), lam20=c(q20_scale, q20_shape),
		     meanlog, sdlog, Q=0, sexactive=1)
    
    saveRDS(list(dat=dat, est=est, swvar=swvar), file=output1)

    #saveRDS(est, file=output1)
    
    print(k)
  }
  
}
