
#source("dt_gen_fct.r")

# generating with all exponentials

#This program generates data for HPV natural history and creates simulated trial data (regular intervals) and cross-sectional data (prevalence study)
#WLOG, we simulate just one HPV type
#LCC - 21Aug2023 allows for precancers to regress to normal

dt_csd_gen.f <- function(n, meanlog, sdlog, cuts, lam01, q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape){
  
  #simulate age of first sex
  fsex <-rlnorm(n,meanlog,sdlog)
  #create natural history (NH); complete markov 
  transitions <- data.frame()
  transitions <- rbind(transitions,t_appear2(fsex, cuts, lam01)) #1st appearance
  transitions <- rbind(transitions,infect(transitions$to, q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape)) #1st infection resolved
  transitions <- rbind(transitions,immune(transitions$to[(nrow(transitions)-n+1):nrow(transitions)])) #natural immunity
  
  #repeat until individuals are too old to be included in studies
  for (i in 0:10){
    transitions <- rbind(transitions,t_appear2(transitions$to[(nrow(transitions)-n+1):nrow(transitions)], cuts, lam01))
    transitions <- rbind(transitions,infect(transitions$to[(nrow(transitions)-n+1):nrow(transitions)], q10_scale, q10_shape, q12_scale, q12_shape, q20_scale, q20_shape))
    transitions <- rbind(transitions,immune(transitions$to[(nrow(transitions)-n+1):nrow(transitions)]))
  }
  sapply(transitions[(nrow(transitions)-n+1):nrow(transitions),],summary) #check of age and state distribution at end of simulated NH
  transitions <- rbind(transitions, data.frame(id=1:n, from=0, to=fsex, from_x=0, to_x=0))
  transitions <- transitions[order(transitions$id,transitions$to),]
  
  print(dim(transitions))
  
  #create survey data
  survey_inclusion <- rep(1, n) #sample(0:1, n, replace = TRUE)
  age_at_survey <-sample(14:60, n, replace = TRUE)
  survey_ids <- seq(n)[survey_inclusion==1]
  
  transitions_select = subset(transitions, id %in% survey_ids)
  dat0=dplyr::left_join(transitions_select,data.frame(id=survey_ids, 
                                                      age=age_at_survey[survey_inclusion==1], 
                                                      fsex=fsex[survey_inclusion==1]))
  dat0$ind=with(dat0, from< age & age <= to)
  
  print(dim(dat0))
  
  
  dat1=subset(dat0, ind==1)
  
  
  surveys=dat1[,c("id", "fsex", "age")]
  surveys$x <- ifelse(dat1$to==dat1$age, dat1_to_x, dat1$from_x)
  surveys$w <- 1
  
  print(dim(surveys))
  
  # --- 
  surveys$hpv <- ifelse(surveys$x==0, 0, 1)
  table(surveys$age, surveys$hpv, useNA = "always")
  dat=subset(surveys, !is.na(surveys$hpv) & (surveys$fsex<surveys$age)) # sextual active 
  print(dim(dat))
  
  return(dat)
}

