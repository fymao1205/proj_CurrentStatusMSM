

P11.f <- function(t, lam10, lam12){
  pweibull(t, lam10[2], lam10[1], 0,0)*pweibull(t, lam12[2], lam12[1], 0, 0)
}

P12.f <- function(t, lam10, lam12, lam20){
  
  # generate gaussian quad matrix
  gau1=gaussleg(max(20, ceiling(t)), 0, t)
  u1=gau1[,1]; w1=gau1[,2]
  
  res1=dweibull(u1, lam12[2], lam12[1])*pweibull(u1, lam10[2], lam10[1], 0, 0)*pweibull(t-u1, lam20[2], lam20[1], 0,0)
  
  res=sum(res1*w1)

  return(res) 
}

P00.f <- function(l, r, cuts_S, lam01_S){ # l > S; r-l <=1
  
  res=exp(-(Hpwc_double(r, cuts_S, lam01_S, 0)-Hpwc_double(l, cuts_S, lam01_S, 0)))
  
  return(res)
}

G.f <- function(r, R, lam10, lam12, lam20){
  
  res1=P11.f(R-r, lam10, lam12)
  res2=P12.f(R-r, lam10, lam12, lam20)
  res=res1+res2
  res=ifelse(is.nan(res)|is.na(res), res1, res)
  return(res)
}

G12.f <- function(r, R, lam10, lam12, lam20){
  
  res1=P11.f(R-r, lam10, lam12)
  res2=P12.f(R-r, lam10, lam12, lam20)
  res2=ifelse(is.nan(res2), 0, res2)
  res=c(res1,res2)
  return(res)
}


# S: age of first sex; double
# R: current age; int
Tk.f <- function(lk, uk, R, cuts_S, lam01_S, lam10, lam12, lam20){
  
  gau1=gaussleg(20, lk, uk)
  u1=gau1[,1]; w1=gau1[,2]
  
  lam01_r=sapply(u1, function(x){hpwc_double(x, cuts_S, lam01_S, 0)})
  P00_r=sapply(u1, function(x){
    exp(-(Hpwc_double(x, cuts_S, lam01_S, 0)-Hpwc_double(lk, cuts_S, lam01_S, 0)))
  })
  f_r=P00_r*lam01_r

  G_r=sapply(u1, function(x){G.f(x, R, lam10, lam12, lam20)})
  
  res=sum(f_r*G_r*w1)

  return(res)
}

Tk12.f <- function(lk, uk, R, cuts_S, lam01_S, lam10, lam12, lam20){
  
  gau1=gaussleg(20, lk, uk)
  u1=gau1[,1]; w1=gau1[,2]
  
  lam01_r=sapply(u1, function(x){hpwc_double(x, cuts_S, lam01_S, 0)})
  P00_r=sapply(u1, function(x){
    #ppwc_double(x, cuts_S, lam01_S, 0, 0)
    exp(-(Hpwc_double(x, cuts_S, lam01_S, 0)-Hpwc_double(lk, cuts_S, lam01_S, 0)))
  })
  f_r=P00_r*lam01_r
  #f_r = ifelse(is.nan(f_r), 0, f_r)
  
  G12_r=lapply(u1, function(x){G12.f(x, R, lam10, lam12, lam20)})
  G_r=do.call("rbind", G12_r)
  #print(dim(G_r))
  
  res=colSums(f_r*G_r*w1)

  return(res)
}

P_Z_g_S.f <- function(S, R, cuts, lam01, lam10, lam12, lam20){
  
  K=ceiling(R)-floor(S)-1
  
  if(any(cuts>S)){
    cuts_S=c(S, cuts[c(cuts)>S])
    lam01_S=c(0, lam01[c(cuts, Inf)>S])
  }else{
    cuts_S=S
    lam01_S=c(0, lam01[length(lam01)])
  }
  
  if(K==0){
    Tvec=Tk.f(S, R, R, cuts_S, lam01_S, lam10, lam12, lam20)
  }else{
    
    lkvec=c(S, floor(S)+1:K)
    ukvec=c(floor(S)+1:K, R)
    
    Tvec=sapply(1:length(lkvec), function(k){
      Tk.f(lkvec[k],ukvec[k], R, cuts_S, lam01_S, lam10, lam12, lam20)
    })
  }
  
  P0vec=c(1)
  for(k in 1:(K+1)){
    p0k=1-sum(P0vec*Tvec[1:k])
    P0vec=c(P0vec, p0k)
  }
  
  P12_g_S_vec=sum(Tvec*P0vec[-(K+2)])
  
  return(list(Tvec=Tvec, P0vec=P0vec, 
              P12_g_S=1-P0vec[(K+2)],
              P12_g_S_alt=P12_g_S_vec))
}

P_Z12_g_S.f <- function(S, R, cuts, lam01, lam10, lam12, lam20){
  
  K=ceiling(R)-floor(S)-1
  
  if(any(cuts>S)){
    cuts_S=c(S, cuts[c(cuts)>S])
    lam01_S=c(0, lam01[c(cuts, Inf)>S])
  }else{
    cuts_S=S
    lam01_S=c(0, lam01[length(lam01)])
  }
  
  if(K==0){
    Tvec=Tk12.f(S, R, R, cuts_S, lam01_S, lam10, lam12, lam20)
    Tmat=matrix(Tvec, ncol=2)
  }else{
    
    lkvec=c(S, floor(S)+1:K)
    ukvec=c(floor(S)+1:K, R)
    
    Tlist=lapply(1:length(lkvec), function(k){
      Tk12.f(lkvec[k],ukvec[k], R, cuts_S, lam01_S, lam10, lam12, lam20)
    })
    Tvec=do.call("rbind", Tlist)
    Tmat=as.matrix(Tvec)
  }
  
  # print(dim(Tmat))
  
  P0vec=c(1)
  for(k in 1:(K+1)){
    p0k=1-sum(P0vec*Tmat[1:k,])
    P0vec=c(P0vec, p0k)
  }
  
  P12_g_S_vec=colSums(Tmat*P0vec[-(K+2)])
  
  return(list(Tvec=Tmat, P0vec=P0vec, 
              Ppos_g_S=1-P0vec[(K+2)],
              P12_g_S=P12_g_S_vec))
}



