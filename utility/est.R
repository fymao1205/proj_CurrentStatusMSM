

# log likelihood given S
logL_g_S.f <- function(dt, cuts, lam01, lam10, lam12, lam20, num.cores){
  
  # res2=apply(dt, 1, function(x){
  #   res1=P_Z_g_S.f(x[5], x[2], cuts, lam01, lam10, lam12, lam20)$P12_g_S
  #   res=ifelse(x[6]==1, res1, 1-res1)
  #   log(res)
  # })
  res2=mclapply(1:nrow(dt), function(x){
    res1=P_Z_g_S.f(dt[x,"fsex"], dt[x,"age"], cuts, lam01, lam10, lam12, lam20)$P12_g_S
    res=ifelse(dt[x,"hpv"]==1, res1, 1-res1)
    log(res)
  }, mc.cores=num.cores)
  
  logL=as.vector(unlist(res2))
  logL=ifelse(is.na(logL), 0, logL)
  logL=ifelse(is.finite(logL), logL, 0)
  
  return(logL)
}

# estimating function given S
est01_g_S.f <- function(ini, dt, cuts, lam10, lam12, lam20, num.cores){
  
  len=length(cuts)+1
  
  if(is.null(ini)){
    ini=c(rep(log(0.01), len+1))#, 0.01)
    #ini=log(0.01)
  }
  
  obj.f<-function(para){
    
    lam01=exp(para[1:len])
    #lam01=c(exp(para), 0.05)
    #lam20=c(exp(-para[len+1]), 1)
    
    res0=logL_g_S.f(dt, cuts, lam01, lam10, lam12, lam20, num.cores)
    res0=ifelse(is.na(res0), 0, res0)
    res=-sum(res0)
    
    print(res)
    if(is.na(res)){
      print(which(is.na(res0)))
    }
    
    print(para)
    return(res)
    
  }
  
  res=optim(par=ini,
            obj.f, method=("L-BFGS-B"), #lower=c(rep(-10, len), 0),#upper=c(rep(5,len), Inf), lower=c(rep(-10, len), 0),
            hessian=TRUE)

  return(res)
  
}


est01_20_g_S.f <- function(ini, dt, cuts, lam10, lam12, num.cores){
  
  len=length(cuts)+1
  
  if(is.null(ini)){
    ini=c(rep(log(0.01), len+1))#, 0.01)
    #ini=log(0.01)
  }
  
  obj.f<-function(para){
    
    lam01=exp(para[1:len])
    #lam01=c(exp(para), 0.05)
    lam20=c(exp(-para[len+1]), 1)
    
    res0=logL_g_S.f(dt, cuts, lam01, lam10, lam12, lam20, num.cores)
    res0=ifelse(is.na(res0), 0, res0)
    res0=ifelse(is.finite(res0), res0, 0)
    res=-sum(res0)
    
    print(res)
    if(is.na(res)){
      print(which(is.na(res0)))
    }
    
    print(para)
    return(res)
    
  }
  
  res=optim(par=ini,
            obj.f, method=("L-BFGS-B"), #lower=c(rep(-10, len), 0),#upper=c(rep(5,len), Inf), lower=c(rep(-10, len), 0),
            hessian=TRUE)
  
  return(res)
  
}


score_01_g_S.f <- function(gradstep=1e-05, par, dt, cuts, lam10, lam12, lam20, num.cores){
  
  para.len=length(par)
  lam01=exp(par[1:len])
  
  logL0=logL_g_S.f(dt, cuts, lam01, lam10, lam12, lam20, num.cores)

  logLp.list <- lapply(1:para.len, function(i){
    print(i)
    par.p=par
    par.p[i]=par[i]+gradstep
    lam01.p=exp(par.p[1:len])
    #loglik_rec20_exp.f(par.p, cuts, lam12, lam10, meanlog, sdlog, Q, q1_vec, dtw)
    logL_g_S.f(dt, cuts, lam01.p, lam10, lam12, lam20, num.cores)
  })
  
  s.list <- lapply(logLp.list, function(x){
    #(x$loglik-logL0$loglik)/gradstep
    (x-logL0)/gradstep
  })
  
  smat=do.call("cbind",s.list)
  
  return(list(score=smat,
              loglik0=logL0,
              loglikp=logLp.list))
  
}

swvar_01_g_S.f <- function(gradstep=1e-05, est.res, dt, cuts, lam10, lam12, lam20, num.cores){
  
  
  len=length(cuts)+1
 
  hess=est.res$hess
  smat=score_01_g_S.f(gradstep, est.res$par, dt, cuts, lam10, lam12, lam20, num.cores)
  
  B = t(smat$score) %*% (smat$score)
  IA = ginv(hess)
  V=IA %*% B %*% IA
  
  ase=sqrt(diag(V))
  
  return(list(B=B, V=V, se.loglam=ase, score=smat))
}


score_01_20_g_S.f <- function(gradstep=1e-05, par, dt, cuts, lam10, lam12, num.cores){
  
  para.len=length(par)
  lam01=exp(par[1:para.len])
  lam20=c(exp(-par[para.len+1]), 1)
  
  logL0=logL_g_S.f(dt, cuts, lam01, lam10, lam12, lam20, num.cores)
  
  logLp.list <- lapply(1:para.len, function(i){
    print(i)
    par.p=par
    par.p[i]=par[i]+gradstep
    lam01.p=exp(par.p[1:para.len])
    lam20.p=c(exp(-par.p[para.len+1]), 1)
    #loglik_rec20_exp.f(par.p, cuts, lam12, lam10, meanlog, sdlog, Q, q1_vec, dtw)
    logL_g_S.f(dt, cuts, lam01.p, lam10, lam12, lam20.p, num.cores)
  })
  
  s.list <- lapply(logLp.list, function(x){
    #(x$loglik-logL0$loglik)/gradstep
    (x-logL0)/gradstep
  })
  
  smat=do.call("cbind",s.list)
  
  return(list(score=smat,
              loglik0=logL0,
              loglikp=logLp.list))
  
}

swvar_01_20_g_S.f <- function(gradstep=1e-05, est.res, dt, cuts, lam10, lam12, num.cores){
  
  
  hess=est.res$hess
  smat=score_01_20_g_S.f(gradstep, est.res$par, dt, cuts, lam10, lam12, num.cores)
  
  B = t(smat$score) %*% (smat$score)
  IA = ginv(hess)
  V=IA %*% B %*% IA
  
  ase=sqrt(diag(V))
  
  return(list(B=B, V=V, se.loglam=ase, score=smat))
}


