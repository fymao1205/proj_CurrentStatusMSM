
# utility functions for data generation

#function to simulate time to HPV appearance
#-- strt_x: current state (mostly 0)
#-- strt: age of first sex

library(eha)

pwc_S.f <- function(S, cuts, lam01){
  
  #lapply(1:length(S), function(i){
  if(any(cuts>S)){
    cuts_S=c(S, cuts[cuts>S])
    lam01_S=c(0, lam01[c(cuts, Inf)>S])
  }else{
    cuts_S=S
    lam01_S=c(0, lam01[length(lam01)])
  }
  res <- list(cuts_S=cuts_S, lam01_S=lam01_S)
  #})
  return(res)
}

t_appear2 <- function(strt, cuts, lam01){
  
  #u_n=runif(length(strt), 0, 1)
  
  t_n <- lapply(1:length(strt), function(i){
    pwc_S.list = pwc_S.f(strt[i], cuts, lam01)
    res <- try(rpch(1, cuts=pwc_S.list$cuts_S, levels=pwc_S.list$lam01_S))
    
    if(isTRUE(class(res)=="try-error")) { NA }else{ res }
    
  })#, mc.cores=2)

  return(data.frame(id=seq(length(strt)),from=strt,to=as.numeric(unlist(t_n)),
                    from_x=rep(0,length(strt)),to_x=rep(1,length(strt))))
}


#function to simulate resolution of HPV infection
infect <- function(strt, q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape){
  clr <- rweibull(length(strt),q10_shape,q10_scale)
  prg <- rweibull(length(strt),q12_shape,q12_scale)
  reg <- rweibull(length(strt),q20_shape,q20_scale) #rgamma(length(strt),q20_shape,scale=q20_scale)
  t_trans <- strt+pmin(clr,prg)
  t_trans <- ifelse(is.na(t_trans), Inf, t_trans)
  x_trans <- ifelse(clr<prg,0,2)
  t_trans2 <- t_trans+reg
  t_trans2 <- ifelse(is.na(t_trans2), Inf, t_trans2)
  #dataset of all transitions
  temp <- rbind(data.frame(id=seq(length(strt)),from=strt,to=t_trans,from_x=rep(1,length(strt)),to_x=x_trans),
                data.frame(id=seq(length(strt)),from=t_trans,to=t_trans2,from_x=rep(2,length(strt)),to_x=rep(0,length(strt)))[which(x_trans==2),])
  #reorder transitions so clearance set occurs last
  temp1 <- subset(temp,to_x==2)
  temp2 <- subset(temp,to_x==0)
  temp2 <- temp2[order(temp2$id),]
  return(rbind(temp1,temp2))
}

#function to simulate natural immunity post-infection
immune <- function(strt){
  t_trans <- strt+d # --- d: length of immune
  return(data.frame(id=seq(length(strt)),from=strt,to=t_trans,from_x=rep(0,length(strt)),to_x=rep(0,length(strt))))
}


dt_csd_gen.f <- function(n, meanlog, sdlog, cuts, lam01, q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape){
  
  #simulate age of first sex
  fsex <-rlnorm(n,meanlog,sdlog)
  #print(paste0("fsex: ",range(fsex)))
  #create natural history (NH); complete markov 
  transitions <- data.frame()
  transitions <- rbind(transitions,t_appear2(fsex, cuts, lam01)) #1st appearance
  #print(paste0("iter0: t_appear2 ", range(transitions$to)))
  transitions <- rbind(transitions,infect(transitions$to, q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape)) #1st infection resolved
  #transitions <- rbind(transitions,immune(transitions$to[(nrow(transitions)-n+1):nrow(transitions)])) #natural immunity
  #print(paste0("iter0: infect ",range(transitions$to)))
  
  #repeat until individuals are too old to be included in studies
  for (i in 0:10){
    transitions <- rbind(transitions,t_appear2(transitions$to[(nrow(transitions)-n+1):nrow(transitions)], cuts, lam01))
    #print(range(transitions$to))
    transitions <- rbind(transitions,infect(transitions$to[(nrow(transitions)-n+1):nrow(transitions)], q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape))
    #transitions <- rbind(transitions,immune(transitions$to[(nrow(transitions)-n+1):nrow(transitions)]))
    #print(range(transitions$to))
  }
  sapply(transitions[(nrow(transitions)-n+1):nrow(transitions),],summary) #check of age and state distribution at end of simulated NH
  transitions <- rbind(transitions, data.frame(id=1:n, from=0, to=fsex, from_x=0, to_x=0))
  transitions <- transitions[order(transitions$id,transitions$to),]
  
  #create survey data
  survey_inclusion <- rep(1, n)#sample(0:1, n, replace = TRUE)
  age_at_survey <-sample(14:60, n, replace = TRUE)
  survey_ids <- seq(n)[survey_inclusion==1]
  
  transitions_select = subset(transitions, id %in% survey_ids)
  dat0=dplyr::left_join(transitions_select,data.frame(id=survey_ids, age=age_at_survey[survey_inclusion==1], 
                                                      fsex=fsex[survey_inclusion==1]))
  dat0$ind=with(dat0, from< age & age <= to)
  dat1=subset(dat0, ind==1)
  
  surveys=dat1[,c("id", "fsex", "age")]
  surveys$x <- ifelse(dat1$to==dat1$age, dat1_to_x, dat1$from_x)
  surveys$w <- 1
  
  # --- 
  surveys$hpv <- ifelse(surveys$x==0, 0, 1)
  table(surveys$age, surveys$hpv, useNA = "always")
  dat=subset(surveys, !is.na(surveys$hpv) & (surveys$fsex<surveys$age)) # sextual active 
  
  return(dat)
}

