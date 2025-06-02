

library(Rcpp)
library(eha)
library(flexsurv)
library(parallel)
library(MASS)
sourceCpp("/home/maof3/HPV_appearance/sim1/utility/commonf.cpp")
source("/home/maof3/HPV_appearance/sim1/utility/dt_gen_fct.R")
#source("/home/maof3/HPV_appearance/sim1/utility/dt_csd_gen.R")
source("/home/maof3/HPV_appearance/sim1/utility/commonf.R")
source("/home/maof3/HPV_appearance/sim1/utility/est.R")

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
#q10_scale <- 0.75
#q10_shape <- 0.74
#weibull parameters for progression to precancer
#q12_scale <- 4.755734#15/2
#q12_shape <- 2.2#1.2
#weibull parameters for regression of precancers
#q20_scale <- 10#2.5#2.5
#q20_shape <- 1#0.4

step1.par=c( log(c(0.75, 0.74)), 0.017522708, 2.2,  log(10), 1)
step1.var=diag( c(0.018, 0.015, 0.0018, 0.31, 1)^2 )


sim.fun <- function(n.sim, output1,seed, b){
  
  set.seed(seed)
  
  for(k in 1:n.sim){
    
    q10_scale <- exp(rnorm(1, mean=log(0.75), sd=0.018))
    q10_shape <- exp(rnorm(1, mean=log(0.74), sd=0.015))
    #weibull parameters for progression to precancer
    # CIN3+ of HPV18/45

    a1 <- rnorm(1, 1.819, 0.111); a2 <- rnorm(1, -0.888, 0.154)
    q12_shape <- 1/exp(a2) #1/rnorm(1, 0.017522708, sd=0.0018)/12 #exp(rnorm(1, mean=log(4.75), sd=0.5))#4.75#15/2
    q12_scale <- exp(a1) #rnorm(1, 2.2, sd=0.31) #exp(rnorm(1, mean=log(2.2), sd=0.5))#2.2#1.2
    #weibull parameters for regression of precancers
    q20_scale <- exp(rnorm(1, mean=log(12.5), sd=1))
    q20_shape <- 1

    dat=dt_csd_gen.f(n, meanlog, sdlog, cuts, lam01, q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape)

   
    # estimation
    #est=est01_g_S.f(ini=rep(0.01, 2), dat, cuts=c(22), c(0.75, 0.74), c(exp(1.819),1/exp(-0.888)),
    #		     c(12.5, 1), num.cores=2)

#    for(b in 1:10){
    set.seed(9999+b)

    bid <- sample(1:nrow(dat), nrow(dat), replace=TRUE)
    dat.b <- dat[bid,]
    q10_scale.b <- exp(rnorm(1, mean=log(0.75), sd=0.018))
    q10_shape.b <- exp(rnorm(1, mean=log(0.74), sd=0.015))
    #weibull parameters for progression to precancer
    #a1.b <- rnorm(1, 2.513, 0.259); a2.b <- rnorm(1, -0.631, 0.192)  
    a1.b <- rnorm(1, 1.819, 0.111); a2.b <- rnorm(1, -0.888, 0.154)
    q12_shape.b <- 1/exp(a2.b) #1/rnorm(1, 0.017522708, sd=0.0018)/12 #exp(rnorm(1, mean=log(4.75), sd=0.5))#4.75#15/2
    q12_scale.b <- exp(a1.b) #rnorm(1, 2.2, sd=0.31) #exp(rnorm(1, mean=log(2.2), sd=0.5))#2.2#1.2
    #weibull parameters for regression of precancers
    q20_scale.b <- exp(rnorm(1, mean=log(12.5), sd=1))##10#2.5#2.5
    q20_shape <- 1
    
    # estimation 
    ini=rep(0.01, 2) #swvar$est$par
    est=est01_g_S.f(ini=ini, dat.b, cuts=c(22), c(q10_scale.b, q10_shape.b), c(q12_scale.b, q12_shape.b),
                    c(q20_scale.b, q20_shape), num.cores=2)
    
    # save the results
    saveRDS(list(dat=dat.b, est=est), file=output1)
    
    #print(k)
  }
  
}
