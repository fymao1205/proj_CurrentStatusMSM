
# likelihood construction and estimation when state at recruitment is known 



# estimating lam01[1]

P11.f <- function(t, lam10, lam12){
  pweibull(t, lam10[2], lam10[1], 0,0)*pweibull(t, lam12[2], lam12[1], 0, 0)
}

P12.f <- function(t, lam10, lam12, lam20){
  
  # generate gaussian quad matrix
  #gau1=gaussleg(20, 0, t)
  gau1=gaussleg(max(20, ceiling(t)), 0, t)
  u1=gau1[,1]; w1=gau1[,2]
  
  #res1=P11.f(u1, lam10, lam12)
  res1=dweibull(u1, lam12[2], lam12[1])*pweibull(u1, lam10[2], lam10[1], 0, 0)*pweibull(t-u1, lam20[2], lam20[1], 0,0)
  
  res=sum(res1*w1)
  #res=sum(res1*w1)
  
  return(res) 
}

P00.f <- function(l, r, cuts_S, lam01_S){ # l > S; r-l <=1
  
  #gau1=gaussleg(20, l, r)
  #u1=gau1[,1]; w1=gau1[,2]
  
  #lam01_r=sapply(u1, function(x){hpwc_double(x, cuts_S, lam01_S, 0)})
  
  res=exp(-(Hpwc_double(r, cuts_S, lam01_S, 0)-Hpwc_double(l, cuts_S, lam01_S, 0)))
  
  return(res)
}

G.f <- function(r, R, lam10, lam12, lam20){
  
  res1=P11.f(R-r, lam10, lam12)
  res2=P12.f(R-r, lam10, lam12, lam20)
  res=res1+res2
  res=ifelse(is.nan(res), res1, res)
  res=ifelse(is.na(res), res1, res)
  #res=c(res1, res2)
  return(res)
}


# S: age of first sex; double
# R: current age; int
Tk.f <- function(lk, uk, R, cuts_S, lam01_S, lam10, lam12, lam20){
  
  gau1=gaussleg(20, lk, uk)
  u1=gau1[,1]; w1=gau1[,2]
  
  #lam01_r=sapply(u1, function(x){hpwc_double(x, c(S,cuts), c(0,lam01), 0)})
  #P00_r=sapply(u1, function(x){ppwc_double(x, c(S,cuts), c(0,lam01), 0, 0)})
  
  # if(any(cuts>S)){
  #   cuts_S=c(S, cuts[c(cuts)>S])
  #   lam01_S=c(0, lam01[c(cuts, Inf)>S])
  # }else{
  #   cuts_S=S
  #   lam01_S=c(0, lam01[length(lam01)])
  # }
  
  # if(any(cuts>lk)){
  #   cuts_S=c(lk, cuts[c(cuts)>lk])
  #   lam01_S=c(0, lam01[c(cuts, Inf)>lk])
  # }else{
  #   cuts_S=lk
  #   lam01_S=c(0, lam01[length(lam01)])
  # }
  
  lam01_r=sapply(u1, function(x){hpwc_double(x, cuts_S, lam01_S, 0)})
  #lam01_r=sapply(u1, function(x){hpch(x, cuts_S, lam01_S, 0)})
  P00_r=sapply(u1, function(x){
    #ppwc_double(x, cuts_S, lam01_S, 0, 0)
    exp(-(Hpwc_double(x, cuts_S, lam01_S, 0)-Hpwc_double(lk, cuts_S, lam01_S, 0)))
    #exp(-(Hpch(x, cuts_S, lam01_S, 0)-Hpch(lk, cuts_S, lam01_S, 0)))
  })
  f_r=P00_r*lam01_r
  #f_r = ifelse(is.nan(f_r), 0, f_r)
  
  G_r=sapply(u1, function(x){G.f(x, R, lam10, lam12, lam20)})
  
  res=sum(f_r*G_r*w1)
  #res=sum(lam01_r*P00_r*G_r*w1)#/ppwc_double(lk, cuts_S, lam01_S, 0, 0)
  
  return(res)
}

T1k.f <- function(lk, uk, R, cuts_S, lam01_S, lam10, lam12, lam20){
  
  gau1=gaussleg(20, lk, uk)
  u1=gau1[,1]; w1=gau1[,2]
  
  #lam01_r=sapply(u1, function(x){hpwc_double(x, c(S,cuts), c(0,lam01), 0)})
  #P00_r=sapply(u1, function(x){ppwc_double(x, c(S,cuts), c(0,lam01), 0, 0)})
  
  # if(any(cuts>S)){
  #   cuts_S=c(S, cuts[c(cuts)>S])
  #   lam01_S=c(0, lam01[c(cuts, Inf)>S])
  # }else{
  #   cuts_S=S
  #   lam01_S=c(0, lam01[length(lam01)])
  # }
  
  # if(any(cuts>lk)){
  #   cuts_S=c(lk, cuts[c(cuts)>lk])
  #   lam01_S=c(0, lam01[c(cuts, Inf)>lk])
  # }else{
  #   cuts_S=lk
  #   lam01_S=c(0, lam01[length(lam01)])
  # }
  
  lam01_r=sapply(u1, function(x){hpwc_double(x, cuts_S, lam01_S, 0)})
  #lam01_r=sapply(u1, function(x){hpch(x, cuts_S, lam01_S, 0)})
  P00_r=sapply(u1, function(x){
    #ppwc_double(x, cuts_S, lam01_S, 0, 0)
    exp(-(Hpwc_double(x, cuts_S, lam01_S, 0)-Hpwc_double(lk, cuts_S, lam01_S, 0)))
    #exp(-(Hpch(x, cuts_S, lam01_S, 0)-Hpch(lk, cuts_S, lam01_S, 0)))
  })
  f_r=P00_r*lam01_r
  #f_r = ifelse(is.nan(f_r), 0, f_r)
  
  #G_r=sapply(u1, function(x){G.f(x, R, lam10, lam12, lam20)})
  G_r=sapply(u1, function(x){P11.f(R-x, lam10, lam12)})
  
  res=sum(f_r*G_r*w1)
  #res=sum(lam01_r*P00_r*G_r*w1)#/ppwc_double(lk, cuts_S, lam01_S, 0, 0)
  
  return(res)
}

T2k.f <- function(lk, uk, R, cuts_S, lam01_S, lam10, lam12, lam20){
  
  gau1=gaussleg(20, lk, uk)
  u1=gau1[,1]; w1=gau1[,2]
  
  #lam01_r=sapply(u1, function(x){hpwc_double(x, c(S,cuts), c(0,lam01), 0)})
  #P00_r=sapply(u1, function(x){ppwc_double(x, c(S,cuts), c(0,lam01), 0, 0)})
  
  # if(any(cuts>S)){
  #   cuts_S=c(S, cuts[c(cuts)>S])
  #   lam01_S=c(0, lam01[c(cuts, Inf)>S])
  # }else{
  #   cuts_S=S
  #   lam01_S=c(0, lam01[length(lam01)])
  # }
  
  # if(any(cuts>lk)){
  #   cuts_S=c(lk, cuts[c(cuts)>lk])
  #   lam01_S=c(0, lam01[c(cuts, Inf)>lk])
  # }else{
  #   cuts_S=lk
  #   lam01_S=c(0, lam01[length(lam01)])
  # }
  
  lam01_r=sapply(u1, function(x){hpwc_double(x, cuts_S, lam01_S, 0)})
  #lam01_r=sapply(u1, function(x){hpch(x, cuts_S, lam01_S, 0)})
  P00_r=sapply(u1, function(x){
    #ppwc_double(x, cuts_S, lam01_S, 0, 0)
    exp(-(Hpwc_double(x, cuts_S, lam01_S, 0)-Hpwc_double(lk, cuts_S, lam01_S, 0)))
    #exp(-(Hpch(x, cuts_S, lam01_S, 0)-Hpch(lk, cuts_S, lam01_S, 0)))
  })
  f_r=P00_r*lam01_r
  #f_r = ifelse(is.nan(f_r), 0, f_r)
  
  #G_r=sapply(u1, function(x){G.f(x, R, lam10, lam12, lam20)})
  G_r=sapply(u1, function(x){P12.f(R-x, lam10, lam12, lam20)})
  
  res=sum(f_r*G_r*w1)
  #res=sum(lam01_r*P00_r*G_r*w1)#/ppwc_double(lk, cuts_S, lam01_S, 0, 0)
  
  return(res)
}

# P(Z(R)=j|S); j=1,2
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
    #Tvec=0
    Tvec=Tk.f(S, R, R, cuts_S, lam01_S, lam10, lam12, lam20)
    T1vec=T1k.f(S, R, R, cuts_S, lam01_S, lam10, lam12, lam20)
    T2vec=T2k.f(S, R, R, cuts_S, lam01_S, lam10, lam12, lam20)
  }else{
    
    lkvec=c(S, floor(S)+1:K)
    ukvec=c(floor(S)+1:K, R)
    
    Tvec=sapply(1:length(lkvec), function(k){
      Tk.f(lkvec[k],ukvec[k], R, cuts_S, lam01_S, lam10, lam12, lam20)
    })
    
    T1vec=sapply(1:length(lkvec), function(k){
      T1k.f(lkvec[k],ukvec[k], R, cuts_S, lam01_S, lam10, lam12, lam20)
    })
    
    T2vec=sapply(1:length(lkvec), function(k){
      T2k.f(lkvec[k],ukvec[k], R, cuts_S, lam01_S, lam10, lam12, lam20)
    })
  }

  P0vec=c(1)
  for(k in 1:(K+1)){
    p0k=1-sum(P0vec*Tvec[1:k])
    P0vec=c(P0vec, p0k)
  }
  
  P12_g_S_vec=sum(Tvec*P0vec[-(K+2)])
    
  P1_g_S_vec=sum(T1vec*P0vec[-(K+2)])
  P2_g_S_vec=sum(T2vec*P0vec[-(K+2)])
  
  P12_g_S_vec=P1_g_S_vec+P2_g_S_vec
  
  return(list(Tvec=Tvec, P0vec=P0vec, 
              P12_g_S=1-P0vec[(K+2)],
              P12_g_S_alt=P12_g_S_vec,
              P1_g_S=P1_g_S_vec,
              P2_g_S=P2_g_S_vec))
}

# matrix: Sp_vec * Rp_vec
Pmat_Z1.f <- function(Rp_vec, cuts, lam01, lam10, lam12, lam20, agof.dist){
  
  Rmax=max(Rp_vec)
  gauS=gaussleg(max(20, ceiling(Rmax)), 0, Rmax)
  Sp_vec=gauS[,1]; wSp=gauS[,2]
  dSp=dgengamma(Sp_vec, agof.dist$meanlog, agof.dist$sdlog, agof.dist$Q)
  
  res.mat=matrix(0, nrow=length(Sp_vec), ncol=length(Rp_vec))
  
  for(r in 1:length(Rp_vec)){
    
    rp=Rp_vec[r]
    #print(rp)
    res=mclapply(Sp_vec[Sp_vec<=rp], function(sp){
      P_Z_g_S.f(sp, rp, cuts, lam01, lam10, lam12, lam20)$P12_g_S
    }, mc.cores=2)
    res.mat[Sp_vec<=rp,r] <- as.vector(unlist(res))
    
  }
  
  P1_Rp=apply(res.mat, 2, function(x){
    sum(x*wSp*dSp)
  })
  
  
  return(list(p1.marg=P1_Rp, 
              p1gs.mat=res.mat,
              dSp=dSp,
              Sp=Sp_vec))
}


Pmat_Z12.f <- function(Rp_vec, cuts, lam01, lam10, lam12, lam20, agof.dist){
  
  Rmax=max(Rp_vec)
  gauS=gaussleg(max(20, ceiling(Rmax)), 0, Rmax)
  Sp_vec=gauS[,1]; wSp=gauS[,2]
  dSp=dgengamma(Sp_vec, agof.dist$meanlog, agof.dist$sdlog, agof.dist$Q)
  
  res2.mat=res1.mat=matrix(0, nrow=length(Sp_vec), ncol=length(Rp_vec))
  
  for(r in 1:length(Rp_vec)){
    
    rp=Rp_vec[r]
    #print(rp)
    res=mclapply(Sp_vec[Sp_vec<=rp], function(sp){
      #P_Z_g_S.f(sp, rp, cuts, lam01, lam10, lam12, lam20)$P12_g_S
      res=P_Z12_g_S.f(sp, rp, cuts, lam01, lam10, lam12, lam20)
      c(res$P1_g_S, res$P2_g_S)
    }, mc.cores=2)
    res0.mat=do.call("rbind", res)
    res1.mat[Sp_vec<=rp,r] <- res0.mat[,1]#as.vector(unlist(res))
    res2.mat[Sp_vec<=rp,r] <- res0.mat[,2]
  }
  
  P1_Rp=apply(res1.mat, 2, function(x){
    sum(x*wSp*dSp)
  })
  
  P2_Rp=apply(res2.mat, 2, function(x){
    sum(x*wSp*dSp)
  })
  
  return(list(p1.marg=P1_Rp, 
              p2.marg=P2_Rp, 
              p1gs.mat=res1.mat,
              p2gs.mat=res2.mat,
              dSp=dSp,
              Sp=Sp_vec))
}

est_lam01_pwc01_wei_rec20_g_j.f <- function(ini, dtw, cuts, lam12, lam10, meanlog, sdlog, Q, sexactive=1){
  
  len=length(cuts)+1
  
  if(is.null(ini)){
    ini=c(rep(log(0.01), len+1))#, 0.01)
    #ini=log(0.01)
  }
  
  age_vec=sort(unique(dtw$age))
  
  if(sexactive==1){
    q1_vec= pgengamma(age_vec, meanlog, sdlog, Q)
  }else{
    q1_vec=1#rep(1, length(age_vec))
  }
  
  
  #create_dt_res=create_dt_list.f(age_vec, cuts)
  #dim(create_dt_res[[2]])
  
  #if(length(create_dt_res[[1]])!=len){stop("data list does not match the pwc!")}
  
  obj.f<-function(para){
    
    lam01=exp(para[1:len])
    #lam01=c(exp(para), 0.05)
    lam20=c(exp(-para[len+1]), 1)
    
    #p1_vec=Pmat_Z1.f(age_vec, cuts, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))$p1.marg
    res=Pmat_Z12.f(age_vec, cuts, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))
    p1_vec=res$p1.marg
    p2_vec=res$p2.marg
    p1_sex=p1_vec/q1_vec
    p2_sex=p2_vec/q1_vec
    p0_sex=1-(p1_sex+p2_sex)
    
    print("done p1_g_L_res ...")
    
    dtp=data.frame(age=age_vec, 
                   p2=p2_sex, 
                   p1=p1_sex, 
                   p0=p0_sex)
    colnames(dtp)=c("age", "p2", "p1", "p0")
    dt0=dplyr::left_join(dtw,dtp)
    dt=dt0[,c("age","x","w")]
    dt$p=ifelse(dt$x==2, dt0$p2, ifelse(dt$x==1, dt0$p1, dt0$p0))
    if(any(dt$p>1) | any(dt$p<0)){ print(paste0("warning: ", dt$p[dt$p>1], dt$p[dt$p<0])) }
    
    dt$p=ifelse(dt$p>1, 1, dt$p)
    dt$p=ifelse(dt$p<0, 0, dt$p)    
    res0=dt$w*log(dt$p);
    #res0=ifelse(dt$p==0, 0, res0)
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
  #res=nlm(f=obj.f, p=ini, hessian=TRUE)
  
  return(list(est=res,
              #datalong=create_dt_res,
              q1_vec=q1_vec))
}

loglik_rec20_exp_g_j.f <- function(par, cuts, lam12, lam10, meanlog, sdlog, Q, q1_vec, dtw){
  
  len=length(cuts)+1
  age_vec=sort(unique(dtw$age))
  
  lam01=exp(par[1:len])
  lam20=c(exp(-par[len+1]),1)
  
  # p1_vec=Pmat_Z12.f(age_vec, cuts, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))$p1.marg
  # 
  # p1_sex=p1_vec/q1_vec
  # p0_sex=1-p1_sex
  # 
  # print("done p1_g_L_res ...")
  # 
  # dtp=data.frame(age=age_vec, 
  #                p1=p1_sex, 
  #                p0=p0_sex)
  # colnames(dtp)=c("age", "p1", "p0")
  # dt0=dplyr::left_join(dtw,dtp)
  # dt=dt0[,c("age","hpv","w")]
  # dt$p=ifelse(dt$hpv==1, dt0$p1, dt0$p0)
  
  res=Pmat_Z12.f(age_vec, cuts, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))
  p1_vec=res$p1.marg
  p2_vec=res$p2.marg
  p1_sex=p1_vec/q1_vec
  p2_sex=p2_vec/q1_vec
  p0_sex=1-(p1_sex+p2_sex)
  
  print("done p1_g_L_res ...")
  
  dtp=data.frame(age=age_vec, 
                 p2=p2_sex, 
                 p1=p1_sex, 
                 p0=p0_sex)
  colnames(dtp)=c("age", "p2", "p1", "p0")
  dt0=dplyr::left_join(dtw,dtp)
  dt=dt0[,c("age","x","w")]
  dt$p=ifelse(dt$x==2, dt0$p2, ifelse(dt$x==1, dt0$p1, dt0$p0))
  if(any(dt$p>1) | any(dt$p<0)){ print(paste0("warning: ", dt$p[dt$p>1], dt$p[dt$p<0])) }
  
  dt$p=ifelse(dt$p>1, 1, dt$p)
  dt$p=ifelse(dt$p<0, 0, dt$p)    
  
  res0=dt$w*log(dt$p);
  
  return(list(loglik=res0, p1=p1_vec))
  
  #return(res)
  
}

score_rec20_exp_g_j.f <- function(gradstep=1e-05, par, cuts, lam12, lam10, meanlog, sdlog, Q, q1_vec, dtw){
  
  para.len=length(par)
  
  logL0=loglik_rec20_exp_g_j.f(par, cuts, lam12, lam10, meanlog, sdlog, Q, q1_vec, dtw)
  
  logLp.list <- lapply(1:para.len, function(i){
    print(i)
    par.p=par
    par.p[i]=par[i]+gradstep
    loglik_rec20_exp_g_j.f(par.p, cuts, lam12, lam10, meanlog, sdlog, Q, q1_vec, dtw)
  })
  
  s.list <- lapply(logLp.list, function(x){
    (x$loglik-logL0$loglik)/gradstep
    #(x-logL0)/gradstep
  })
  
  smat=do.call("cbind",s.list)
  
  return(list(score=smat,
              loglik0=logL0,
              loglikp=logLp.list))
  
}

swvar_rec20_exp_g_j.f <- function(gradstep=1e-05, est.res, dtw, cuts, lam12, lam10, meanlog, sdlog, Q, sexactive=1){
  
  
  len=length(cuts)+1
  age_vec=sort(unique(dtw$age))
  
  if(sexactive==1){
    q1_vec= pgengamma(age_vec, meanlog, sdlog, Q)
  }else{
    q1_vec=1#rep(1, length(age_vec))
  }
  
  #create_dt_res=create_dt_list.f(age_vec, cuts)
  #dim(create_dt_res[[2]])
  
  hess=est.res$hess
  smat=score_rec20_exp_g_j.f(gradstep, est.res$est$par, cuts, lam12, lam10, meanlog, sdlog, Q, q1_vec, dtw)
  
  B = t(smat$score) %*% (smat$score/ifelse(dtw$w==0, 1, dtw$w))
  IA = ginv(hess)
  V=IA %*% B %*% IA
  
  ase=sqrt(diag(V))
  
  # d_p1_mat
  p1_0 <- smat$loglik0$p1
  p1_p.list <- lapply(smat$loglikp, function(x){x$p1})
  d_p1.list <- lapply(p1_p.list, function(x){
    (x-p1_0)/gradstep
  })
  
  d_p1.mat=do.call("cbind", d_p1.list)
  
  var_p1 <- apply(d_p1.mat, 1, function(x){
    t(x) %*% V %*% (x)
  })
  
  # se_p1 <- apply(d_p1.mat, 1, function(x){
  #   t(x) %*% ase
  # })
  
  se_p1 = sqrt(as.vector(unlist(var_p1)))
  
  return(list(B=B, V=V, se.loglam=ase, se.lam=ase*exp(est.res$est), se.p1=se_p1, score=smat))
}



# log likelihood given S
loglik_rec20_exp_g_j_S.f <- function(dt, cuts, lam01, lam10, lam12, lam20, num.cores){
  
  # res2=apply(dt, 1, function(x){
  #   res1=P_Z_g_S.f(x[5], x[2], cuts, lam01, lam10, lam12, lam20)$P12_g_S
  #   res=ifelse(x[6]==1, res1, 1-res1)
  #   log(res)
  # })
  res2=mclapply(1:nrow(dt), function(x){
    res1=P_Z12_g_S.f(dt[x,"fsex"], dt[x,"age"], cuts, lam01, lam10, lam12, lam20)
    res=ifelse(dt[x,"x"]==1, res1$P1_g_S, ifelse(dt[x,"x"]==2, res1$P2_g_S, 1-res1$P1_g_S-res1$P2_g_S))
    log(res)
  }, mc.cores=num.cores)
  
  logL=as.vector(unlist(res2))
  logL=ifelse(is.na(logL), 0, logL)
  logL=ifelse(is.finite(logL), logL, 0)
  
  return(logL)
}

est_lam01_pwc01_wei_rec20_g_j_S.f <- function(ini, dt, cuts, lam12, lam10, num.scores){
  
  len=length(cuts)+1
  
  if(is.null(ini)){
    ini=c(rep(log(0.01), len+1))#, 0.01)
    #ini=log(0.01)
  }
  
  obj.f<-function(para){
    
    lam01=exp(para[1:len])
    #lam01=c(exp(para), 0.05)
    lam20=c(exp(-para[len+1]), 1)
    
    res0=loglik_rec20_exp_g_j_S.f(dt, cuts, lam01, lam10, lam12, lam20, num.cores)
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
  #res=nlm(f=obj.f, p=ini, hessian=TRUE)
  
  return(list(est=res))
}


score_rec20_exp_g_j_S.f <- function(gradstep=1e-05, par, cuts, lam12, lam10, dt, num.cores){
  
  par.len=length(par)
  
  len=length(cuts)+1
  
  lam01=exp(par[1:len])
  lam20=c(exp(-par[len+1]), 1)
  
  logL0=loglik_rec20_exp_g_j_S.f(dt, cuts, lam01, lam10, lam12, lam20, num.cores)
  
  logLp.list <- lapply(1:par.len, function(i){
    print(i)
    par.p=par
    par.p[i]=par[i]+gradstep
    lam01.p=exp(par.p[1:len])
    lam20.p=c(exp(-par.p[len+1]), 1)
    loglik_rec20_exp_g_j_S.f(dt, cuts, lam01.p, lam10, lam12, lam20.p, num.cores)
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


swvar_rec20_exp_g_j_S.f <- function(gradstep=1e-05, est.res, cuts, lam12, lam10, dt, num.cores){
  
  
  len=length(cuts)+1
  #create_dt_res=create_dt_list.f(age_vec, cuts)
  #dim(create_dt_res[[2]])
  
  hess=est.res$est$hess
  smat=score_rec20_exp_g_j_S.f(gradstep, est.res$est$par, cuts, lam12, lam10, dt, num.cores)
  
  B = t(smat$score) %*% (smat$score)
  IA = ginv(hess)
  V=IA %*% B %*% IA
  
  ase=sqrt(diag(V))
  
  # d_p1_mat
  #p1_0 <- smat$loglik0$p1
  #p1_p.list <- lapply(smat$loglikp, function(x){x$p1})
  #d_p1.list <- lapply(p1_p.list, function(x){
  #  (x-p1_0)/gradstep
  #})
  
  #d_p1.mat=do.call("cbind", d_p1.list)
  
  #var_p1 <- apply(d_p1.mat, 1, function(x){
  #  t(x) %*% V %*% (x)
  #})
  
  # se_p1 <- apply(d_p1.mat, 1, function(x){
  #   t(x) %*% ase
  # })
  
  #se_p1 = sqrt(as.vector(unlist(var_p1)))
  
  return(list(B=B, V=V, se.loglam=ase, se.lam=ase*exp(est.res$est$par), #se.p1=se_p1, 
              score=smat))
}


# test=est_lam01_pwc01_wei_rec20.f(rep(0.01,3), dtw, cuts, lam12, lam10, meanlog, sdlog, Q, 1)
# 
# 
# 
# Rp_vec= unique(dtw$age); Rmax=max(Rp_vec)
# gauS=gaussleg(max(20, ceiling(Rmax)), 0, Rmax)
# Sp_vec=gauS[,1]
#   
# 
# P_Z_g_S.f(13.5, 15, 22, lam01, lam10, lam12, lam20)
# P_Z12_g_S.f(13.5, 15, 22, lam01, lam10, lam12, lam20)
# 
# P_Z_g_S.f(13.5, 35, lam01, lam10, lam12, lam20)
# 
# P_Z_g_S.f(13.5, 14, cuts, lam01, lam10, lam12, lam20)
# 
# 
#Pmat_Z1.f(15:17, 22, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))$p1.marg
#test=Pmat_Z12.f(15:17, 22, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))
# Pmat_Z1.f(15:17, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))$p1gs.mat
# Pmat_Z1.f(15:17, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))$dSp
# Pmat_Z1.f(15:17, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))$Sp
# 
# 
# test=Pmat_Z1.f(Rp_vec, lam01, lam10, lam12, lam20, list(meanlog=meanlog, sdlog=sdlog, Q=Q))

