

library(Rcpp)
library(eha)
library(flexsurv)
library(parallel)
sourceCpp("./utility/commonf.cpp")
source("./utility/dt_gen_fct.r")
source("./utility/est_rec20.r")

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

#lam10=c(q10_scale, q10_shape)
#lam12=c(q12_scale, q12_shape)
#lam20=c(q20_scale, q20_shape)

sim.fun <- function(n.sim, output1,seed){
  
  set.seed(seed)

  for(k in 1:n.sim){
    
    #set.seed(k)
    # data generation
    fsex <-rlnorm(n,meanlog,sdlog)
  
  dt.list <- lapply(1:n, function(i){
    h_gen.f(i, fsex[i], lam01, cuts, lam10=c(q10_scale, q10_shape), lam12=c(q12_scale, q12_shape), lam20=c(q20_scale, q20_shape))
  })
  dt=do.call("rbind", dt.list)
  dt=rbind(dt, data.frame(id=1:n, from=0, to=fsex, from_x=0, to_x=0))
  dt=dt[order(dt$id,dt$to),]
  
  # check 1
  if(any(dt$to < dt$from)){stop("history generation wrong")}
  
  #create survey data
  survey_inclusion <- sample(0:1, n, replace = TRUE)
  #survey_inclusion <- rbinom(n, 1, z)
  age_at_survey <-sample(14:60, n, replace = TRUE)#runif(n,14,60)
  survey_ids <- seq(n)[survey_inclusion==1]
  
  dt_select = subset(dt, id %in% survey_ids)
  dat0=dplyr::left_join(dt_select,data.frame(id=survey_ids, age=age_at_survey[survey_inclusion==1]))
  dat0$ind=with(dat0, from< age & age <= to)
  dat1=subset(dat0, ind==1)
  
  # check 2
  if(nrow(dat1) != sum(survey_inclusion)){stop("selection wrong")}
  
  surveys=dat1[,c("id", "age")]
  surveys$x <- ifelse(dat1$to==dat1$age, dat1_to_x, dat1$from_x)
  surveys$w <- 1#1/z[survey_inclusion]
  
  # --- 
  surveys$fsex=fsex[survey_inclusion==1]
  surveys$hpv <- ifelse(surveys$x==0, 0, 1)
  table(surveys$age, surveys$hpv, useNA = "always")
  
  dat=subset(surveys, !is.na(hpv) & fsex<=age ) # sextual active 
  any(dat$fsex>dat$age & dat$hpv==1)
  
  dtw<-aggregate(w~age+hpv,dat,sum)
    
    # estimation
    ini=c(rep(log(0.01),3))
    est=est_lam01_pwc01_wei_rec20.f(ini, dtw, cuts=c(22), lam12=c(q12_scale, q12_shape), lam10=c(q10_scale, q10_shape), meanlog, sdlog, Q=0, sexactive=1)

    saveRDS(c(est$est$est), file=output1)
    
    print(k)
  }
  
}
