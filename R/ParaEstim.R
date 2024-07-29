Estimate=function(
    para,  # True parameter
    n,     # sample size
    m,     # No. of failure times
    L,     # Lower specifiaction limit
    U,     # Upper specifiaction limit
    T,     # Target value
    GM,     # gamma value
    MC.size  =11000,  # Size of MCMC sample   
    MC.burn  =1000,   # Burn-in sample
    q.val=-0.5, # genarlized entropy loss parameter
    c.val=0.5,  # linex loss parameter
    Sim.no=500, # No. pof replications
    randomseed=c(2021,6,30), #seed nomber
    direc   # directory saving the results
    ){
  cat("\014")
  library(Rcpp)
  library(pracma)
  library(roptim)
  library(VGAM)
  library(coda)
  hessian=numDeriv::hessian  # to use hessian in numDerive package
  
  setwd(direc)

  options(width=100,length=200,max.print = 10000,digits = 3, scipen=999)

  ## Censroing schemes
  Scheme=matrix(nr=3,nc=m,0)
  Scheme[1,]=c(rep(0,m-1),n-m)
  Scheme[2,]=c(n-m,rep(0,m-1))
  Scheme[3,]=c(rep(1,n-m),rep(0,2*m-n))
  
  nu=para[1]
  xi=para[2]
  ini.para=para
  
  spmk=spmk_fun(para,L,U,T,GM)
  para=c(para,spmk)
  for(i1 in 1:3){
    R=Scheme[i1,]  
    set.seed(randomseed)
    it=1
    while (it<=Sim.no){
      if(it==1){
        MLE=MPS=EM=NULL
        MH_LK_SEL=MH_LK_GEL=MH_LK_LIN=NULL
        MH_PS_SEL=MH_PS_GEL=MH_PS_LIN=NULL
        ARS_LK_SEL=ARS_LK_GEL=ARS_LK_LIN=NULL
        ARS_PS_SEL=ARS_PS_GEL=ARS_PS_LIN=NULL
        TK_LK_SEL=TK_LK_GEL=TK_LK_LIN=NULL
        TK_PS_SEL=TK_PS_GEL=TK_PS_LIN=NULL
        start.time=date()
      }  
      
      X=Gen_Data(nu,xi,R)
      if(is.character(X)) next
      
      #################################################################
      #------------- MLE using NR method
      ############s#####################################################
      mle.res=Estim(para,X,R,L,U,T,GM,ini.para,"LK",c(0.51,0.1),2*para[1:2],"Nelder-Mead")
      # Check convergence
      mle.res
      if(is.character(mle.res)) next
    
      #################################################################
      #------------- MPS using NR method
      ############s#####################################################
      mps.res=suppressWarnings(try(Estim(para,X,R,L,U,T,GM,ini.para,"PS",c(0.51,0),3*para[1:2],"L-BFGS-B"),silent=TRUE))
      
      # Check convergence
      if(is.character(mps.res)) next
      
      #################################################################
      #------------- MLE using EM algorithm
      ############s#####################################################
      em.res=suppressWarnings(try(EM_Alg(para,X,R,L,U,T,GM,ini.para,upper=100,MaxIter=100,tol=0.0001, verbose=0),silent=TRUE))
      
      # Check convergence
      if(is.character(em.res)) next
      
      mle.est=as.numeric(mle.res$Par)
      mle.bias=as.numeric(mle.res$Bias)
      mle.ese=sqrt(as.numeric(mle.res$Var))
      mle.aci=mle.res$ACI
      
      mps.est=as.numeric(mps.res$Par)
      mps.bias=as.numeric(mps.res$Bias)
      mps.ese=sqrt(as.numeric(mps.res$Var))
      mps.aci=mps.res$ACI
      
      em.est=as.numeric(em.res$Par)
      em.bias=as.numeric(em.res$Bias)
      em.ese=sqrt(as.numeric(em.res$Var))
      em.aci=em.res$ACI
      
      MLE =rbind(MLE,c(mle.est[1],mle.bias[1],mle.aci[1,],
                       mle.est[2],mle.bias[2],mle.aci[2,],
                       mle.est[3],mle.bias[3],mle.aci[3,]))
      MPS =rbind(MPS,c(mps.est[1],mps.bias[1],mps.aci[1,],
                       mps.est[2],mps.bias[2],mps.aci[2,],
                       mps.est[3],mps.bias[3],mps.aci[3,]))
      EM =rbind(EM,c(em.est[1],em.bias[1],em.aci[1,],
                     em.est[2],em.bias[2],em.aci[2,],
                     em.est[3],em.bias[3],em.aci[3,]))
      
      #################################################################
      #------------- Bayes estimate using TK method based on LK function
      ############s#####################################################
      tk.lk.res=TK(para,X,R,L,U,T,GM,q.val,c.val,mle.est[1:2],"LK")
      TK_LK_SEL =rbind(TK_LK_SEL,c(tk.lk.res$est_sel[1],tk.lk.res$bias_sel[1],NA,NA,NA,NA,
                                   tk.lk.res$est_sel[2],tk.lk.res$bias_sel[2],NA,NA,NA,NA,
                                   tk.lk.res$est_sel[3],tk.lk.res$bias_sel[3],NA,NA,NA,NA))
      TK_LK_GEL =rbind(TK_LK_GEL,c(tk.lk.res$est_gel[1],tk.lk.res$bias_gel[1],NA,NA,NA,NA,
                                   tk.lk.res$est_gel[2],tk.lk.res$bias_gel[2],NA,NA,NA,NA,
                                   tk.lk.res$est_gel[3],tk.lk.res$bias_gel[3],NA,NA,NA,NA))
      TK_LK_LIN =rbind(TK_LK_LIN,c(tk.lk.res$est_lin[1],tk.lk.res$bias_lin[1],NA,NA,NA,NA,
                                   tk.lk.res$est_lin[2],tk.lk.res$bias_lin[2],NA,NA,NA,NA,
                                   tk.lk.res$est_lin[3],tk.lk.res$bias_lin[3],NA,NA,NA,NA))
      
      #################################################################
      #------------- Bayes estimate using TK method based on PS function
      ############s#####################################################
      tk.ps.res=TK(para,X,R,L,U,T,GM,q.val,c.val,mle.est[1:2],"PS")
      TK_PS_SEL =rbind(TK_PS_SEL,c(tk.ps.res$est_sel[1],tk.ps.res$bias_sel[1],NA,NA,NA,NA,
                                   tk.ps.res$est_sel[2],tk.ps.res$bias_sel[2],NA,NA,NA,NA,
                                   tk.ps.res$est_sel[3],tk.ps.res$bias_sel[3],NA,NA,NA,NA))
      TK_PS_GEL =rbind(TK_PS_GEL,c(tk.ps.res$est_gel[1],tk.ps.res$bias_gel[1],NA,NA,NA,NA,
                                   tk.ps.res$est_gel[2],tk.ps.res$bias_gel[2],NA,NA,NA,NA,
                                   tk.ps.res$est_gel[3],tk.ps.res$bias_gel[3],NA,NA,NA,NA))
      TK_PS_LIN =rbind(TK_PS_LIN,c(tk.ps.res$est_lin[1],tk.ps.res$bias_lin[1],NA,NA,NA,NA,
                                   tk.ps.res$est_lin[2],tk.ps.res$bias_lin[2],NA,NA,NA,NA,
                                   tk.ps.res$est_lin[3],tk.ps.res$bias_lin[3],NA,NA,NA,NA))
      
      #####################################################################
      #------------- Bayes estimate using MH method based on LK function
      ############s########################################################
      cat("\n Generating MH samples based on likelihood function for iteration ",it,"\n")
      MH_LK=MH_sample(para,X,R,L,U,T,GM,q.val,c.val,mle.est[1:2],"LK",MC.size,MC.burn,mle.ese[1:2],3,F,T)
      if(is.infinite(MH_LK$est_sel[1])) next
      if(is.infinite(MH_LK$est_sel[2])) next
      if(is.infinite(MH_LK$est_sel[3])) next
      if(is.infinite(MH_LK$est_gel[1])) next
      if(is.infinite(MH_LK$est_gel[2])) next
      if(is.infinite(MH_LK$est_gel[3])) next
      if(is.infinite(MH_LK$est_lin[1])) next
      if(is.infinite(MH_LK$est_lin[2])) next
      if(is.infinite(MH_LK$est_lin[3])) next
      
      MH_LK_SEL =rbind(MH_LK_SEL ,c(MH_LK$est_sel[1],MH_LK$bias_sel[1],MH_LK$HPD[1,],
                                    MH_LK$est_sel[2],MH_LK$bias_sel[2],MH_LK$HPD[2,],
                                    MH_LK$est_sel[3],MH_LK$bias_sel[3],MH_LK$HPD[3,]))
      MH_LK_GEL =rbind(MH_LK_GEL ,c(MH_LK$est_gel[1],MH_LK$bias_gel[1],NA,NA,NA,NA,
                                    MH_LK$est_gel[2],MH_LK$bias_gel[2],NA,NA,NA,NA,
                                    MH_LK$est_gel[3],MH_LK$bias_gel[3],NA,NA,NA,NA))
      MH_LK_LIN =rbind(MH_LK_LIN ,c(MH_LK$est_lin[1],MH_LK$bias_lin[1],NA,NA,NA,NA,
                                    MH_LK$est_lin[2],MH_LK$bias_lin[2],NA,NA,NA,NA,
                                    MH_LK$est_lin[3],MH_LK$bias_lin[3],NA,NA,NA,NA))
      
      
      #####################################################################
      #------------- Bayes estimate using MH method based on PS function
      ############s########################################################
      cat("\n Generating MH samples based on product of spacing function for iteration ",it,"\n");
      MH_PS=MH_sample(para,X,R,L,U,T,GM,q.val,c.val,mle.est[1:2],"PS",MC.size,MC.burn,mle.ese[1:2],3,F,T)
      if(is.infinite(MH_PS$est_sel[1])) next
      if(is.infinite(MH_PS$est_sel[2])) next
      if(is.infinite(MH_PS$est_sel[3])) next
      if(is.infinite(MH_PS$est_gel[1])) next
      if(is.infinite(MH_PS$est_gel[2])) next
      if(is.infinite(MH_PS$est_gel[3])) next
      if(is.infinite(MH_PS$est_lin[1])) next
      if(is.infinite(MH_PS$est_lin[2])) next
      if(is.infinite(MH_PS$est_lin[3])) next
      
      
      MH_PS_SEL =rbind(MH_PS_SEL ,c(MH_PS$est_sel[1],MH_PS$bias_sel[1],MH_PS$HPD[1,],
                                    MH_PS$est_sel[2],MH_PS$bias_sel[2],MH_PS$HPD[2,],
                                    MH_PS$est_sel[3],MH_PS$bias_sel[3],MH_PS$HPD[3,]))
      MH_PS_GEL =rbind(MH_PS_GEL ,c(MH_PS$est_gel[1],MH_PS$bias_gel[1],NA,NA,NA,NA,
                                    MH_PS$est_gel[2],MH_PS$bias_gel[2],NA,NA,NA,NA,
                                    MH_PS$est_gel[3],MH_PS$bias_gel[3],NA,NA,NA,NA))
      MH_PS_LIN =rbind(MH_PS_LIN ,c(MH_PS$est_lin[1],MH_PS$bias_lin[1],NA,NA,NA,NA,
                                    MH_PS$est_lin[2],MH_PS$bias_lin[2],NA,NA,NA,NA,
                                    MH_PS$est_lin[3],MH_PS$bias_lin[3],NA,NA,NA,NA))
      
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
      
      sel_mle_tk=round(c(mean(TK_LK_SEL[,1],na.rm=TRUE),
                         mean(TK_LK_SEL[,2],na.rm=TRUE),
                         MSE(nu,TK_LK_SEL[,1],"SEL",0),
                         mean(TK_LK_SEL[,7],na.rm=TRUE),
                         mean(TK_LK_SEL[,8],na.rm=TRUE),
                         MSE(xi,TK_LK_SEL[,7],"SEL",0),
                         mean(TK_LK_SEL[,13],na.rm=TRUE),
                         mean(TK_LK_SEL[,14],na.rm=TRUE),
                         MSE(spmk,TK_LK_SEL[,13],"SEL",0)),3)
      
      gel_mle_tk=round(c(mean(TK_LK_GEL[,1],na.rm=TRUE),
                         mean(TK_LK_GEL[,2],na.rm=TRUE),
                         MSE(nu,TK_LK_GEL[,1],"GEL",q.val),
                         mean(TK_LK_GEL[,7],na.rm=TRUE),
                         mean(TK_LK_GEL[,8],na.rm=TRUE),
                         MSE(xi,TK_LK_GEL[,7],"GEL",q.val),
                         mean(TK_LK_GEL[,13],na.rm=TRUE),
                         mean(TK_LK_GEL[,14],na.rm=TRUE),
                         MSE(spmk,TK_LK_GEL[,13],"GEL",q.val)),3)
      
      lin_mle_tk=round(c(mean(TK_LK_LIN[,1],na.rm=TRUE),
                         mean(TK_LK_LIN[,2],na.rm=TRUE),
                         MSE(nu,TK_LK_LIN[,1],"LINEX",c.val),
                         mean(TK_LK_LIN[,7],na.rm=TRUE),
                         mean(TK_LK_LIN[,8],na.rm=TRUE),
                         MSE(xi,TK_LK_LIN[,7],"LINEX",c.val),
                         mean(TK_LK_LIN[,13],na.rm=TRUE),
                         mean(TK_LK_LIN[,14],na.rm=TRUE),
                         MSE(spmk,TK_LK_LIN[,13],"LINEX",c.val)),3)
      
      sel_mps_tk=round(c(mean(TK_PS_SEL[,1],na.rm=TRUE),
                         mean(TK_PS_SEL[,2],na.rm=TRUE),
                         MSE(nu,TK_PS_SEL[,1],"SEL",0),
                         mean(TK_PS_SEL[,7],na.rm=TRUE),
                         mean(TK_PS_SEL[,8],na.rm=TRUE),
                         MSE(xi,TK_PS_SEL[,7],"SEL",0),
                         mean(TK_PS_SEL[,13],na.rm=TRUE),
                         mean(TK_PS_SEL[,14],na.rm=TRUE),
                         MSE(spmk,TK_PS_SEL[,13],"SEL",0)),3)
      
      gel_mps_tk=round(c(mean(TK_PS_GEL[,1],na.rm=TRUE),
                         mean(TK_PS_GEL[,2],na.rm=TRUE),
                         MSE(nu,TK_PS_GEL[,1],"GEL",q.val),
                         mean(TK_PS_GEL[,7],na.rm=TRUE),
                         mean(TK_PS_GEL[,8],na.rm=TRUE),
                         MSE(xi,TK_PS_GEL[,7],"GEL",q.val),
                         mean(TK_PS_GEL[,13],na.rm=TRUE),
                         mean(TK_PS_GEL[,14],na.rm=TRUE),
                         MSE(spmk,TK_PS_GEL[,13],"GEL",q.val)),3)
      
      lin_mps_tk=round(c(mean(TK_PS_LIN[,1],na.rm=TRUE),
                         mean(TK_PS_LIN[,2],na.rm=TRUE),
                         MSE(nu,TK_PS_LIN[,1],"LINEX",c.val),
                         mean(TK_PS_LIN[,7],na.rm=TRUE),
                         mean(TK_PS_LIN[,8],na.rm=TRUE),
                         MSE(xi,TK_PS_LIN[,7],"LINEX",c.val),
                         mean(TK_PS_LIN[,13],na.rm=TRUE),
                         mean(TK_PS_LIN[,14],na.rm=TRUE),
                         MSE(spmk,TK_PS_LIN[,13],"LINEX",c.val)),3)
      
      sel_mle_mh=round(c(mean(MH_LK_SEL[,1],na.rm=TRUE),
                         mean(MH_LK_SEL[,2],na.rm=TRUE),
                         MSE(nu,MH_LK_SEL[,1],"SEL",0),
                         mean(MH_LK_SEL[,7],na.rm=TRUE),
                         mean(MH_LK_SEL[,8],na.rm=TRUE),
                         MSE(xi,MH_LK_SEL[,7],"SEL",0),
                         mean(MH_LK_SEL[,13],na.rm=TRUE),
                         mean(MH_LK_SEL[,14],na.rm=TRUE),
                         MSE(spmk,MH_LK_SEL[,13],"SEL",0)),3)
      sel_mps_mh=round(c(mean(MH_PS_SEL[,1],na.rm=TRUE),
                         mean(MH_PS_SEL[,2],na.rm=TRUE),
                         MSE(nu,MH_PS_SEL[,1],"SEL",0),
                         mean(MH_PS_SEL[,7],na.rm=TRUE),
                         mean(MH_PS_SEL[,8],na.rm=TRUE),
                         MSE(xi,MH_PS_SEL[,7],"SEL",0),
                         mean(MH_PS_SEL[,13],na.rm=TRUE),
                         mean(MH_PS_SEL[,14],na.rm=TRUE),
                         MSE(spmk,MH_PS_SEL[,13],"SEL",0)),3)
      
      gel_mle_mh=round(c(mean(MH_LK_GEL[,1],na.rm=TRUE),
                         mean(MH_LK_GEL[,2],na.rm=TRUE),
                         MSE(nu,MH_LK_GEL[,1],"GEL",q.val),
                         mean(MH_LK_GEL[,7],na.rm=TRUE),
                         mean(MH_LK_GEL[,8],na.rm=TRUE),
                         MSE(xi,MH_LK_GEL[,7],"GEL",q.val),
                         mean(MH_LK_GEL[,13],na.rm=TRUE),
                         mean(MH_LK_GEL[,14],na.rm=TRUE),
                         MSE(spmk,MH_LK_GEL[,13],"GEL",q.val)),3)
      
      gel_mps_mh=round(c(mean(MH_PS_GEL[,1],na.rm=TRUE),
                         mean(MH_PS_GEL[,2],na.rm=TRUE),
                         MSE(nu,MH_PS_GEL[,1],"GEL",q.val),
                         mean(MH_PS_GEL[,7],na.rm=TRUE),
                         mean(MH_PS_GEL[,8],na.rm=TRUE),
                         MSE(xi,MH_PS_GEL[,7],"GEL",q.val),
                         mean(MH_PS_GEL[,13],na.rm=TRUE),
                         mean(MH_PS_GEL[,14],na.rm=TRUE),
                         MSE(spmk,MH_PS_GEL[,13],"GEL",q.val)),3)
      lin_mle_mh=round(c(mean(MH_LK_LIN[,1],na.rm=TRUE),
                         mean(MH_LK_LIN[,2],na.rm=TRUE),
                         MSE(nu,MH_LK_LIN[,1],"LINEX",c.val),
                         mean(MH_LK_LIN[,7],na.rm=TRUE),
                         mean(MH_LK_LIN[,8],na.rm=TRUE),
                         MSE(xi,MH_LK_LIN[,7],"LINEX",c.val),
                         mean(MH_LK_LIN[,13],na.rm=TRUE),
                         mean(MH_LK_LIN[,14],na.rm=TRUE),
                         MSE(spmk,MH_LK_LIN[,13],"LINEX",c.val)),3)
      
      lin_mps_mh=round(c(mean(MH_PS_LIN[,1],na.rm=TRUE),
                         mean(MH_PS_LIN[,2],na.rm=TRUE),
                         MSE(nu,MH_PS_LIN[,1],"LINEX",c.val),
                         mean(MH_PS_LIN[,7],na.rm=TRUE),
                         mean(MH_PS_LIN[,8],na.rm=TRUE),
                         MSE(xi,MH_PS_LIN[,7],"LINEX",c.val),
                         mean(MH_PS_LIN[,13],na.rm=TRUE),
                         mean(MH_PS_LIN[,14],na.rm=TRUE),
                         MSE(spmk,MH_PS_LIN[,13],"LINEX",c.val)),3)
      
      
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
      colnames(Estimate)=c("nu.Est","nu.Bias","nu.MSE","xi.Est","xi.Bias","xi.MSE","Spmk.Est","Spmk.Bias","Spmk.MSE")
      rownames(Estimate)=c("MLE","MPS","EM", "TK.LK.SEL","TK.LK.GEL","TK.LK.LINEX","TK.PS.SEL","TK.PS.GEL","TK.PS.LINEX","MH.LK.SEL","MH.LK.GEL","MH.LK.LINEX","MH.PS.SEL","MH.PS.GEL","MH.PS.LINEX")
      colnames(ACI)=c("nu.L","nu.U","nu.Len","nu.CP","xi.L","xi.U","xi.Len","xi.CP","Spmk.L","Spmk.U","Spmk.Len","Spmk.CP")
      rownames(ACI)=c("MLE","MPS","EM","MH.LK","MH.PS")
      
      
      cat("\n======= Simulation Results from iteration 1 to iteration ",it,"============\n")
      
      print(round(Estimate,3))
      print(round(ACI,3))
      
      if(it<Sim.no){it=it+1; next} 
  
      save.image(paste(n,"_",m,"_",i1,sep=""))
      sink(paste("Estimate_ND_n",n,".txt",sep=""),append = TRUE)
      cat("====================================================","\n")
      cat("True Value of Epsilon      : ",para[1],"\n")
      cat("True Value of Eta          : ",para[2],"\n")
      cat("Sample events (n)          : ",m,"\n")
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
      print(round(Estimate,3))
      cat("HPD intervals\n")
      print(round(ACI,3))
      sink()
      break
    }
  }
}
