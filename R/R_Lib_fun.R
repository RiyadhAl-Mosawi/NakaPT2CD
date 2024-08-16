## Generatind sample from Nakagami PT2C
Gen_Data=function(para,R){
  m=length(R)
  w=u=v=t=vv=RR=numeric(m)
  w=runif(m,0,1)
  for(i in 1:m){
    RR[i]=0
    for(j in (m-i+1):m) RR[i]=RR[i]+R[j]
    v[i]=w[i]^(1/(i+RR[i]))
  }
  for(i in 1:m){
    vv[i]=1
    for(j in (m-i+1):m)
      vv[i]=vv[i]*v[j]
    u[i]=1-vv[i]
    t[i]=qnaka(u[i],scale=para[2],shape=para[1])
    #t[i]=inv.fun(u[i],para)
  }
  return(sort(t))
}

#-- Metrepolis-Hasting and adaptive Metrepolis algorithm
RAM_Metropolis <- function(type,     # "LK" for likelihood function and "PS" for product space function
                           para,   # initial values of parameters
                           S,      # initial values of var-cov matrix
                           n_iter, # MCMC sample size
                           n_burnin, # MCMC sample burnin size
                           adapt = FALSE,  # adapted or ordinary MH
                           showProgressBar=interactive()) {
  
  ## =======================================================
  ## (Adaptive) Metropolis Sampler
  ## Implementation of the RAM (robust adaptive Metropolis)
  ## sampler of Vihola, M. (2011) Robust adaptive Metropolis algorithm with
  ## coerced acceptance rate. Statistics and Computing.
  post=function(x,Type=LK){
    ll=1
    if(Type=="LK"){
      ll=prod(dnaka(X,scale=x[2],shape=x[1])*(1-pnaka(X,scale=x[2],shape=x[1]))^R)*sqrt(x[1]*trigamma(x[1])-1)/x[2]
    }
    if(Type=="PS"){
      ll=pnaka(X[1],scale=x[2],shape=x[1])*(1-pnaka(X[m],scale=x[2],shape=x[1]))*prod(pnaka(X[2:m],scale=x[2],shape=x[1])-pnaka(X[1:(m-1)],scale=x[2],shape=x[1]))
      ll=ll*prod((1-pnaka(X,scale=x[2],shape=x[1]))^R)**sqrt(x[1]*trigamma(x[1])-1)/x[2]
    }
    return(ll)
  }
  
  p <- length(para)
  chain <- matrix(NA, n_iter, p)
  accept <- numeric(n_iter)
  chain[1, ] <- para
  posterior <- log(post(para,Type=type))
  
  #cat('  generate', n_iter, 'samples \n')
  if(showProgressBar){
    pb <- txtProgressBar(min=0, max=n_iter, style=3)
  }
  
  for (i in 2:n_iter){
    if(showProgressBar) {
      setTxtProgressBar(pb, i)
    }
    
    u <- rnorm(p)
    para_prop <- chain[i - 1, ] + S %*% u
    if (!any(para_prop<=0)) {
      acceptance_prob <- min(1, exp(log(post(para_prop,Type=type)) - log(post(chain[i - 1, ],Type=type))))
      if(is.nan(acceptance_prob)) acceptance_prob=0
      if (runif(1) < acceptance_prob) {
        accept[i] <- 1
        chain[i, ] <- para_prop
      }else{
        chain[i, ] <- chain[i - 1, ]
      }
    } else {
      chain[i, ] <- chain[i - 1, ]
      acceptance_prob <- 0
    }
    if(adapt & i <= n_burnin) {
      S <- ramcmc::adapt_S(S, u, acceptance_prob, i - 1)
    }
  }
  if(showProgressBar){
    close(pb)                             # close progress bar
  }
  list(chain = chain[(n_burnin + 1):n_iter, ], S = S,
       acceptance_rate = sum(accept[(n_burnin + 1):n_iter]) / (n_iter - n_burnin))
}

MCMC_MH=function(True_par,X,R,l,u,t,gm,q,c,para,type,MC_size,MC_burn,se,fact,verbose=false,display_progress=true,Adapt=FALSE){
  MH_sel=MH_gel=MH_lin=matrix(nr=MC_size-MC_burn,nc=3)
  HPD=matrix(nr=3,nc=4,0)
  est_sel=est_gel=est_lin=numeric(3)
  bias_sel=bias_gel=bias_lin=numeric(3)
  out=RAM_Metropolis(type,para,S=diag(se^2),n_iter=MC_size, n_burnin=MC_burn, adapt = TRUE,showProgressBar=interactive()) 

  MH_sel[,1:2]=out$chain[,1:2]
  MH_sel[,3]  =apply(as.matrix(1:(MC_size-MC_burn)),1,function(w) spmk_fun(c(MH_sel[w,1],MH_sel[w,2]),L,U,T,GM))
  MH_gel[,1] =apply(as.matrix(1:(MC_size-MC_burn)),1,function(w) MH_sel[w,1]^(-q));
  MH_gel[,2] =apply(as.matrix(1:(MC_size-MC_burn)),1,function(w) MH_sel[w,2]^(-q));
  MH_gel[,3] =apply(as.matrix(1:(MC_size-MC_burn)),1,function(w) MH_sel[w,3]^(-q));
  MH_lin[,1] =apply(as.matrix(1:(MC_size-MC_burn)),1,function(w) exp(-c*MH_sel[w,1]));
  MH_lin[,2] =apply(as.matrix(1:(MC_size-MC_burn)),1,function(w) exp(-c*MH_sel[w,2]));
  MH_lin[,3] =apply(as.matrix(1:(MC_size-MC_burn)),1,function(w) exp(-c*MH_sel[w,3]));
  est_sel=c(mean(MH_sel[,1]),mean(MH_sel[,2]),mean(MH_sel[,3]))
  est_gel=c(mean(MH_gel[,1])^(-1/q),mean(MH_gel[,2])^(-1/q),mean(MH_gel[,3])^(-1/q))
  est_lin=c((-1/c)*log(mean(MH_lin[,1])),(-1/c)*log(mean(MH_lin[,2])),(-1/c)*log(mean(MH_lin[,3])))
  bias_sel=est_sel-True_par
  bias_gel=est_gel-True_par
  bias_lin=est_lin-True_par
  HPD[1,1:3]=c(HPDinterval(as.mcmc(MH_sel[,1]))[1],HPDinterval(as.mcmc(MH_sel[,1]))[2],
               HPDinterval(as.mcmc(MH_sel[,1]))[2]-HPDinterval(as.mcmc(MH_sel[,1]))[1])
  if((HPD[1,1]<=True_par[1])&(True_par[1]<=HPD[1,2])) HPD[1,4]=1 else HPD[1,4]=0
  HPD[2,1:3]=c(HPDinterval(as.mcmc(MH_sel[,2]))[1],HPDinterval(as.mcmc(MH_sel[,2]))[2],
               HPDinterval(as.mcmc(MH_sel[,2]))[2]-HPDinterval(as.mcmc(MH_sel[,2]))[1])
  if((HPD[2,1]<=True_par[2])&(True_par[2]<=HPD[2,2])) HPD[2,4]=1 else HPD[2,4]=0
  HPD[3,1:3]=c(HPDinterval(as.mcmc(MH_sel[,3]))[1],HPDinterval(as.mcmc(MH_sel[,3]))[2],
               HPDinterval(as.mcmc(MH_sel[,3]))[2]-HPDinterval(as.mcmc(MH_sel[,3]))[1])
  if((HPD[3,1]<=True_par[3]) & (True_par[3]<=HPD[3,2])) HPD[3,4]=1 else HPD[3,4]=0
  
  return(list(MH_sel=MH_sel,MH_gel=MH_gel,MH_lin=MH_lin,
              est_sel=est_sel,bias_sel=bias_sel,
              est_gel=est_gel,bias_gel=bias_gel,
              est_lin=est_lin,bias_lin=bias_lin,
              HPD=HPD))
}  

Print_Aggr=function(Iter_No=i1, SaveResults=FALSE){
  #---------------------------------------------------------------------
  # Print the results from the first iteration to the current iteration
  #--------------------------------------------------------------------
  Estimate=NULL
  ACI=NULL
  
  mle=mps=em=sel=gel=lin=NULL
  mle=c(mean(MLE[,1],na.rm=TRUE),
        mean(MLE[,2],na.rm=TRUE),
        mean(MLE[,2]^2,na.rm=TRUE),
        mean(MLE[,7],na.rm=TRUE),
        mean(MLE[,8],na.rm=TRUE),
        mean(MLE[,8]^2,na.rm=TRUE),
        mean(MLE[,13],na.rm=TRUE),
        mean(MLE[,14],na.rm=TRUE),
        mean(MLE[,14]^2,na.rm=TRUE))
  
  mps=c(mean(MPS[,1],na.rm=TRUE),
        mean(MPS[,2],na.rm=TRUE),
        mean(MPS[,2]^2,na.rm=TRUE),
        mean(MPS[,7],na.rm=TRUE),
        mean(MPS[,8],na.rm=TRUE),
        mean(MPS[,8]^2,na.rm=TRUE),
        mean(MPS[,13],na.rm=TRUE),
        mean(MPS[,14],na.rm=TRUE),
        mean(MPS[,14]^2,na.rm=TRUE))
  
  em=c(mean(EM[,1],na.rm=TRUE),
       mean(EM[,2],na.rm=TRUE),
       mean(EM[,2]^2,na.rm=TRUE),
       mean(EM[,7],na.rm=TRUE),
       mean(EM[,8],na.rm=TRUE),
       mean(EM[,8]^2,na.rm=TRUE),
       mean(EM[,13],na.rm=TRUE),
       mean(EM[,14],na.rm=TRUE),
       mean(EM[,14]^2,na.rm=TRUE))
  
  sel_mle_tk=c(mean(TK_LK_SEL[,1],na.rm=TRUE),
               mean(TK_LK_SEL[,2],na.rm=TRUE),
               MSE(nu,TK_LK_SEL[,1],"SEL",0),
               mean(TK_LK_SEL[,7],na.rm=TRUE),
               mean(TK_LK_SEL[,8],na.rm=TRUE),
               MSE(xi,TK_LK_SEL[,7],"SEL",0),
               mean(TK_LK_SEL[,13],na.rm=TRUE),
               mean(TK_LK_SEL[,14],na.rm=TRUE),
               MSE(spmk,TK_LK_SEL[,13],"SEL",0))
  
  gel_mle_tk=c(mean(TK_LK_GEL[,1],na.rm=TRUE),
               mean(TK_LK_GEL[,2],na.rm=TRUE),
               MSE(nu,TK_LK_GEL[,1],"GEL",q.val),
               mean(TK_LK_GEL[,7],na.rm=TRUE),
               mean(TK_LK_GEL[,8],na.rm=TRUE),
               MSE(xi,TK_LK_GEL[,7],"GEL",q.val),
               mean(TK_LK_GEL[,13],na.rm=TRUE),
               mean(TK_LK_GEL[,14],na.rm=TRUE),
               MSE(spmk,TK_LK_GEL[,13],"GEL",q.val))
  
  lin_mle_tk=c(mean(TK_LK_LIN[,1],na.rm=TRUE),
               mean(TK_LK_LIN[,2],na.rm=TRUE),
               MSE(nu,TK_LK_LIN[,1],"LINEX",c.val),
               mean(TK_LK_LIN[,7],na.rm=TRUE),
               mean(TK_LK_LIN[,8],na.rm=TRUE),
               MSE(xi,TK_LK_LIN[,7],"LINEX",c.val),
               mean(TK_LK_LIN[,13],na.rm=TRUE),
               mean(TK_LK_LIN[,14],na.rm=TRUE),
               MSE(spmk,TK_LK_LIN[,13],"LINEX",c.val))
  
  sel_mps_tk=c(mean(TK_PS_SEL[,1],na.rm=TRUE),
               mean(TK_PS_SEL[,2],na.rm=TRUE),
               MSE(nu,TK_PS_SEL[,1],"SEL",0),
               mean(TK_PS_SEL[,7],na.rm=TRUE),
               mean(TK_PS_SEL[,8],na.rm=TRUE),
               MSE(xi,TK_PS_SEL[,7],"SEL",0),
               mean(TK_PS_SEL[,13],na.rm=TRUE),
               mean(TK_PS_SEL[,14],na.rm=TRUE),
               MSE(spmk,TK_PS_SEL[,13],"SEL",0))
  
  gel_mps_tk=c(mean(TK_PS_GEL[,1],na.rm=TRUE),
               mean(TK_PS_GEL[,2],na.rm=TRUE),
               MSE(nu,TK_PS_GEL[,1],"GEL",q.val),
               mean(TK_PS_GEL[,7],na.rm=TRUE),
               mean(TK_PS_GEL[,8],na.rm=TRUE),
               MSE(xi,TK_PS_GEL[,7],"GEL",q.val),
               mean(TK_PS_GEL[,13],na.rm=TRUE),
               mean(TK_PS_GEL[,14],na.rm=TRUE),
               MSE(spmk,TK_PS_GEL[,13],"GEL",q.val))
  
  lin_mps_tk=c(mean(TK_PS_LIN[,1],na.rm=TRUE),
               mean(TK_PS_LIN[,2],na.rm=TRUE),
               MSE(nu,TK_PS_LIN[,1],"LINEX",c.val),
               mean(TK_PS_LIN[,7],na.rm=TRUE),
               mean(TK_PS_LIN[,8],na.rm=TRUE),
               MSE(xi,TK_PS_LIN[,7],"LINEX",c.val),
               mean(TK_PS_LIN[,13],na.rm=TRUE),
               mean(TK_PS_LIN[,14],na.rm=TRUE),
               MSE(spmk,TK_PS_LIN[,13],"LINEX",c.val))
  
  sel_mle_mh=c(mean(MH_LK_SEL[,1],na.rm=TRUE),
               mean(MH_LK_SEL[,2],na.rm=TRUE),
               MSE(nu,MH_LK_SEL[,1],"SEL",0),
               mean(MH_LK_SEL[,7],na.rm=TRUE),
               mean(MH_LK_SEL[,8],na.rm=TRUE),
               MSE(xi,MH_LK_SEL[,7],"SEL",0),
               mean(MH_LK_SEL[,13],na.rm=TRUE),
               mean(MH_LK_SEL[,14],na.rm=TRUE),
               MSE(spmk,MH_LK_SEL[,13],"SEL",0))
  sel_mps_mh=c(mean(MH_PS_SEL[,1],na.rm=TRUE),
               mean(MH_PS_SEL[,2],na.rm=TRUE),
               MSE(nu,MH_PS_SEL[,1],"SEL",0),
               mean(MH_PS_SEL[,7],na.rm=TRUE),
               mean(MH_PS_SEL[,8],na.rm=TRUE),
               MSE(xi,MH_PS_SEL[,7],"SEL",0),
               mean(MH_PS_SEL[,13],na.rm=TRUE),
               mean(MH_PS_SEL[,14],na.rm=TRUE),
               MSE(spmk,MH_PS_SEL[,13],"SEL",0))
  
  gel_mle_mh=c(mean(MH_LK_GEL[,1],na.rm=TRUE),
               mean(MH_LK_GEL[,2],na.rm=TRUE),
               MSE(nu,MH_LK_GEL[,1],"GEL",q.val),
               mean(MH_LK_GEL[,7],na.rm=TRUE),
               mean(MH_LK_GEL[,8],na.rm=TRUE),
               MSE(xi,MH_LK_GEL[,7],"GEL",q.val),
               mean(MH_LK_GEL[,13],na.rm=TRUE),
               mean(MH_LK_GEL[,14],na.rm=TRUE),
               MSE(spmk,MH_LK_GEL[,13],"GEL",q.val))
  
  gel_mps_mh=c(mean(MH_PS_GEL[,1],na.rm=TRUE),
               mean(MH_PS_GEL[,2],na.rm=TRUE),
               MSE(nu,MH_PS_GEL[,1],"GEL",q.val),
               mean(MH_PS_GEL[,7],na.rm=TRUE),
               mean(MH_PS_GEL[,8],na.rm=TRUE),
               MSE(xi,MH_PS_GEL[,7],"GEL",q.val),
               mean(MH_PS_GEL[,13],na.rm=TRUE),
               mean(MH_PS_GEL[,14],na.rm=TRUE),
               MSE(spmk,MH_PS_GEL[,13],"GEL",q.val))
  lin_mle_mh=c(mean(MH_LK_LIN[,1],na.rm=TRUE),
               mean(MH_LK_LIN[,2],na.rm=TRUE),
               MSE(nu,MH_LK_LIN[,1],"LINEX",c.val),
               mean(MH_LK_LIN[,7],na.rm=TRUE),
               mean(MH_LK_LIN[,8],na.rm=TRUE),
               MSE(xi,MH_LK_LIN[,7],"LINEX",c.val),
               mean(MH_LK_LIN[,13],na.rm=TRUE),
               mean(MH_LK_LIN[,14],na.rm=TRUE),
               MSE(spmk,MH_LK_LIN[,13],"LINEX",c.val))
  
  lin_mps_mh=c(mean(MH_PS_LIN[,1],na.rm=TRUE),
               mean(MH_PS_LIN[,2],na.rm=TRUE),
               MSE(nu,MH_PS_LIN[,1],"LINEX",c.val),
               mean(MH_PS_LIN[,7],na.rm=TRUE),
               mean(MH_PS_LIN[,8],na.rm=TRUE),
               MSE(xi,MH_PS_LIN[,7],"LINEX",c.val),
               mean(MH_PS_LIN[,13],na.rm=TRUE),
               mean(MH_PS_LIN[,14],na.rm=TRUE),
               MSE(spmk,MH_PS_LIN[,13],"LINEX",c.val))
  
  
  Estimate=data.frame(rbind(mle,mps,em,sel_mle_tk,sel_mps_tk,gel_mle_tk,gel_mps_tk,lin_mle_tk,lin_mps_tk,sel_mle_mh,sel_mps_mh,gel_mle_mh,gel_mps_mh,lin_mle_mh,lin_mps_mh))
  
  aci_mle=aci_mps=aci_em=hpd_mle_mh=hpd_mps_mh=hpd_mle_ars=hpd_mps_ars=NULL
  
  aci_mle=c(max(0,mean(MLE[,3],na.rm=TRUE)),
            mean(MLE[,4],na.rm=TRUE),
            mean(MLE[,5],na.rm=TRUE),
            mean(MLE[,6],na.rm=TRUE),
            max(0,mean(MLE[,9],na.rm=TRUE)),
            mean(MLE[,10],na.rm=TRUE),
            mean(MLE[,11],na.rm=TRUE),
            mean(MLE[,12],na.rm=TRUE),
            max(0,mean(MLE[,15],na.rm=TRUE)),
            mean(MLE[,16],na.rm=TRUE),
            mean(MLE[,17],na.rm=TRUE),
            mean(MLE[,18],na.rm=TRUE))
  aci_mps=c(max(0,mean(MPS[,3],na.rm=TRUE)),
            mean(MPS[,4],na.rm=TRUE),
            mean(MPS[,5],na.rm=TRUE),
            mean(MPS[,6],na.rm=TRUE),
            max(0,mean(MPS[,9],na.rm=TRUE)),
            mean(MPS[,10],na.rm=TRUE),
            mean(MPS[,11],na.rm=TRUE),
            mean(MPS[,12],na.rm=TRUE),
            max(0,mean(MPS[,15],na.rm=TRUE)),
            mean(MPS[,16],na.rm=TRUE),
            mean(MPS[,17],na.rm=TRUE),
            mean(MPS[,18],na.rm=TRUE))
  
  aci_em=c(max(0,mean(EM[,3],na.rm=TRUE)),
           mean(EM[,4],na.rm=TRUE),
           mean(EM[,5],na.rm=TRUE),
           mean(EM[,6],na.rm=TRUE),
           max(0,mean(EM[,9],na.rm=TRUE)),
           mean(EM[,10],na.rm=TRUE),
           mean(EM[,11],na.rm=TRUE),
           mean(EM[,12],na.rm=TRUE),
           max(0,mean(EM[,15],na.rm=TRUE)),
           mean(EM[,16],na.rm=TRUE),
           mean(EM[,17],na.rm=TRUE),
           mean(EM[,18],na.rm=TRUE))
  
  hpd_mle_mh=c(max(0,mean(MH_LK_SEL[,3],na.rm=TRUE)),
               mean(MH_LK_SEL[,4],na.rm=TRUE),
               mean(MH_LK_SEL[,5],na.rm=TRUE),
               mean(MH_LK_SEL[,6],na.rm=TRUE),
               max(0,mean(MH_LK_SEL[,9],na.rm=TRUE)),
               mean(MH_LK_SEL[,10],na.rm=TRUE),
               mean(MH_LK_SEL[,11],na.rm=TRUE),
               mean(MH_LK_SEL[,12],na.rm=TRUE),
               max(0,mean(MH_LK_SEL[,15],na.rm=TRUE)),
               mean(MH_LK_SEL[,16],na.rm=TRUE),
               mean(MH_LK_SEL[,17],na.rm=TRUE),
               mean(MH_LK_SEL[,18],na.rm=TRUE))
  hpd_mps_mh=c(max(0,mean(MH_PS_SEL[,3],na.rm=TRUE)),
               mean(MH_PS_SEL[,4],na.rm=TRUE),
               mean(MH_PS_SEL[,5],na.rm=TRUE),
               mean(MH_PS_SEL[,6],na.rm=TRUE),
               max(0,mean(MH_PS_SEL[,9],na.rm=TRUE)),
               mean(MH_PS_SEL[,10],na.rm=TRUE),
               mean(MH_PS_SEL[,11],na.rm=TRUE),
               mean(MH_PS_SEL[,12],na.rm=TRUE),
               max(0,mean(MH_PS_SEL[,15],na.rm=TRUE)),
               mean(MH_PS_SEL[,16],na.rm=TRUE),
               mean(MH_PS_SEL[,17],na.rm=TRUE),
               mean(MH_PS_SEL[,18],na.rm=TRUE))
  
  
  ACI=data.frame(rbind(aci_mle,aci_mps,aci_em,hpd_mle_mh,hpd_mps_mh))
  colnames(Estimate)=c("nu.Est ","nu.Bias","nu.MSE","xi.Est ","xi.Bias","xi.MSE","Spmk.Est","Spmk.Bias","Spmk.MSE")
  rownames(Estimate)=c("MLE","MPS","EM", "TK.LK.SEL","TK.LK.GEL","TK.LK.LINEX","TK.PS.SEL","TK.PS.GEL","TK.PS.LINEX","MH.LK.SEL","MH.LK.GEL","MH.LK.LINEX","MH.PS.SEL","MH.PS.GEL","MH.PS.LINEX")
  colnames(ACI)=c("nu.L","nu.U","nu.Len","nu.CP","xi.L","xi.U","xi.Len","xi.CP","Spmk.L","Spmk.U","Spmk.Len","Spmk.CP")
  rownames(ACI)=c("MLE","MPS","EM","MH.LK","MH.PS")
  
  
  cat("\n======= Simulation Results from iteration 1 to iteration ",it,"============\n")
  
  print(Estimate,digits =4)
  print(ACI,digits =4)
  
  if(SaveResults==TRUE){
    sink(paste("SimulationResults/Estim_ND_n",n,"_m",m,".txt",sep=""),append = TRUE)
    cat("====================================================","\n")
    cat("True Value of nu           : ",para[1],"\n")
    cat("True Value of xi           : ",para[2],"\n")
    cat("Sample events (n)          : ",n,"\n")
    cat("Number of events (m)       : ",m,"\n")
    cat("No. of iterations          : ",Sim.no,"\n")
    cat("Probability of censrong (I): ",i1,"\n")
    cat("The value of R             : ",R,"\n")
    cat("MCMC sample size           : ",MC.size,"\n")
    cat("MCMC burn-in sample size   : ",MC.burn,"\n")
    cat("GEL loss function para     : ",q.val,"\n")
    cat("LINEX loss function para   : ",c.val,"\n")
    cat("============================================","\n\n")
    cat("Estiamted values\n")
    print(Estimate,digits =4)
    cat("HPD intervals\n")
    print(ACI,digits =4)
    sink()
  }
}



