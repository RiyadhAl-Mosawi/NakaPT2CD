#include <cmath>  // std::pow
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <roptim.h>
// [[Rcpp::depends(roptim)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
// [[Rcpp::depends(RcppProgress)]]

// Library functions of the paper:
// Classical and Bayesian Inference of S'pmk for  Nakagami Distribution Based On PT2C sample

#include<iostream>
#include<algorithm>
#include<algorithm>
#include "gamma_fun.h"
#include "nakagami.h"

using namespace Rcpp;
using namespace std;
using namespace Numer;
using namespace arma;
using namespace roptim;

// These codes for the progress bar inside the loop
#include <progress.hpp>
#include <progress_bar.hpp>

#include "LibFun.h"

  
// [[Rcpp::export]]
Rcpp::List Estim(arma::vec True_Par,
                 arma::vec X, arma::vec R, 
                 double l, double u, double t, double gm,
                 arma::vec para,std::string type, 
                 arma::vec lw, arma::vec up,
                 std::string method) {
  
  LogObjective rb(X,R,type);
  Roptim<LogObjective> opt(method);
  opt.control.trace = 0;
  opt.set_hessian(true);
  opt.set_lower(lw);
  opt.set_upper(up); 
  opt.minimize(rb, para);
  double par_nu=opt.par()[0];
  double par_xi=opt.par()[1];
  double par_spmk=spmk_fun({par_nu,par_xi},l,u,t,gm);
  Rcpp::NumericVector par={par_nu,par_xi,par_spmk};
  Rcpp::NumericMatrix hess  = as<NumericMatrix>(wrap(opt.hessian()));
  Rcpp::NumericMatrix CI(3,4);
  arma::mat v=arma::inv(opt.hessian());
  double var_nu=v(0,0);
  double var_xi=v(1,1);
  Rcpp::NumericVector gr=spmk_grad({par_nu,par_xi},l,u,t,gm);
  double var_spmk=gr(0)*gr(0)*v(0,0)+2*gr(0)*gr(1)*v(0,1)+gr(1)*gr(1)*v(1,1);
  Rcpp::NumericVector var={var_nu,var_xi,var_spmk};
  Rcpp::NumericVector bias={par(0)-True_Par(0),par(1)-True_Par(1),par(2)-True_Par(2)};
   
  CI(0,0)=std::max(0.0,par(0)-1.96*std::sqrt(var(0)));
  CI(0,1)=par(0)+1.96*std::sqrt(var(0));
  CI(0,2)=par(0)+1.96*std::sqrt(var(0))-std::max(0.0,par(0)-1.96*std::sqrt(var(0)));
  CI(0,3)=0;
  if(CI(0,0)<=True_Par(0) && True_Par(0)<=CI(0,1)) CI(0,3)=1.0; 
   
  CI(1,0)=std::max(0.0,par(1)-1.96*std::sqrt(var(1)));
  CI(1,1)=par(1)+1.96*std::sqrt(var(1));
  CI(1,2)=par(1)+1.96*std::sqrt(var(1))-std::max(0.0,par(1)-1.96*std::sqrt(var(1)));
  CI(1,3)=0;
  if(CI(1,0)<=True_Par(1) && True_Par(1)<=CI(1,1)) CI(1,3)=1.0;
   
  CI(2,0)=std::max(0.0,par(2)-1.96*std::sqrt(var(2)));
  CI(2,1)=par(2)+1.96*std::sqrt(var(2));
  CI(2,2)=par(2)+1.96*std::sqrt(var(2))-std::max(0.0,par(2)-1.96*std::sqrt(var(2)));
  CI(2,3)=0;
  if(CI(2,0)<=True_Par(2) && True_Par(2)<=CI(2,1)) CI(2,3)=1.0;
  return Rcpp::List::create(
   Rcpp::Named("Par") = par,
   Rcpp::Named("Bias") = bias,
   Rcpp::Named("Var") = var,
   Rcpp::Named("ACI") = CI,
   Rcpp::Named("Inform") = opt.hessian(),
   Rcpp::Named("Value") = opt.value(),
   Rcpp::Named("Conv") = opt.convergence());
};  

// [[Rcpp::export]]
Rcpp::List EM_Alg(Rcpp::NumericVector True_Par, 
              Rcpp::NumericVector X,
              Rcpp::NumericVector R,
              double l, double u, double t, double gm, 
              Rcpp::NumericVector para, 
              double upper=100, int MaxIter=100, 
              double tol=0.0001, int verbose=0){
  int m=R.size();
  int n=m+sum(R);
  double nu = para[0];
  double xi = para[1];
  
  double dif=10000;
  
  double err_est0, err_est1, err_est2;
  int err_code0, err_code1, err_code2;
  
  Rcpp::NumericVector mle={nu,xi};
  double value=0;
  int it=1;
  double val=0;
  Rcpp::NumericVector init={0.51,2*para[1]};
  while(dif>tol&&it<=MaxIter){
    E0 f0(nu, xi);
    E1 f1(nu, xi);
    E2 f2(nu, xi);
    
    arma::vec W0(m);
    arma::vec W1(m);
    arma::vec W2(m);
    arma::vec Z1(m);
    arma::vec Z2(m);
    
    for(int i=0;i<m;i++){
      W0(i) = integrate(f0, X(i), upper, err_est0, err_code0);
      W1(i) = integrate(f1, X(i), upper, err_est1, err_code1);
      W2(i) = integrate(f2, X(i), upper, err_est2, err_code2);
      Z1(i)=W1(i)/W0(i);     // E(log(X)|X>c)
      Z2(i)=W2(i)/W0(i);    // E(X^2|X>c)
    };
    Rcpp::NumericVector ZZ1=Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(Z1));
    Rcpp::NumericVector ZZ2=Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(Z2));
    xi=0.0;
    for(int i=0;i<m;i++){ 
      xi+=(X(i)*X(i)+R(i)*ZZ2(i))/n;
    };
    Rcpp::Environment base("package:stats"); 
    Rcpp::Function em = base["uniroot"];
    Rcpp::List res = em(Rcpp::_["f"]     = Rcpp::InternalFunction(&cscore),
                        Rcpp::_["interval"]= init,
                        Rcpp::_["xi"]=xi,
                        Rcpp::_["X"]=X,
                        Rcpp::_["R"]=R,
                        Rcpp::_["Z"]=ZZ1);
    nu=res(0);
    dif=abs(nu-mle[0])+abs(xi-mle[1]);
    it+=1;
    mle(0)=nu;
    mle(1)=xi;
    if(verbose>0){
      printf("%f %f %f\n", mle[0], mle[1], dif); 
    };
  }
  Rcpp::NumericVector Par={mle[0],mle[1],spmk_fun(mle,l,u,t,gm)};
  Rcpp::NumericMatrix hess  = inform(mle,X,R,0,upper);
  arma::mat v=arma::inv(as<arma::mat>(wrap(hess)));
  double var_nu=v(0,0);
  double var_xi=v(1,1);
  Rcpp::NumericVector gr=spmk_grad(mle,l,u,t,gm);
  double var_spmk=gr(0)*gr(0)*v(0,0)+2*gr(0)*gr(1)*v(0,1)+gr(1)*gr(1)*v(1,1);
  Rcpp::NumericVector var={var_nu,var_xi,var_spmk};
  Rcpp::NumericVector bias={Par(0)-True_Par(0),Par(1)-True_Par(1),Par(2)-True_Par(2)};
  
  Rcpp::NumericMatrix CI(3,4);
  CI(0,0)=std::max(0.0,Par(0)-1.96*std::sqrt(var(0)));
  CI(0,1)=Par(0)+1.96*std::sqrt(var(0));
  CI(0,2)=Par(0)+1.96*std::sqrt(var(0))-std::max(0.0,Par(0)-1.96*std::sqrt(var(0)));
  CI(0,3)=0;
  if(CI(0,0)<=True_Par(0) && True_Par(0)<=CI(0,1)) CI(0,3)=1.0; 
  
  CI(1,0)=std::max(0.0,Par(1)-1.96*std::sqrt(var(1)));
  CI(1,1)=Par(1)+1.96*std::sqrt(var(1));
  CI(1,2)=Par(1)+1.96*std::sqrt(var(1))-std::max(0.0,Par(1)-1.96*std::sqrt(var(1)));
  CI(1,3)=0;
  if(CI(1,0)<=True_Par(1) && True_Par(1)<=CI(1,1)) CI(1,3)=1.0;
  
  CI(2,0)=std::max(0.0,Par(2)-1.96*std::sqrt(var(2)));
  CI(2,1)=Par(2)+1.96*std::sqrt(var(2));
  CI(2,2)=Par(2)+1.96*std::sqrt(var(2))-std::max(0.0,Par(2)-1.96*std::sqrt(var(2)));
  CI(2,3)=0;
  if(CI(2,0)<=True_Par(2) && True_Par(2)<=CI(2,1)) CI(2,3)=1.0;
  
  return Rcpp::List::create(
    Rcpp::Named("Par") = Par,
    Rcpp::Named("Bias") = bias,
    Rcpp::Named("Var") = var,
    Rcpp::Named("ACI") = CI,
    Rcpp::Named("value") = loglike(mle,X,R),
    Rcpp::Named("inform") = inform(mle,X,R,0,upper));
};  


// This is the main function of TK method
// [[Rcpp::export]]
Rcpp::List TK(arma::vec True_para,
              arma::vec X, 
              arma::vec R, 
              double L, double U, double T,double GM,
              double q, double c,
              arma::vec para, 
              std::string type){
  
  Rcpp::NumericVector est_sel(3),bias_sel(3);
  Rcpp::NumericVector est_gel(3),bias_gel(3);
  Rcpp::NumericVector est_lin(3),bias_lin(3);
  
  int m=R.size();
  int N=sum(R)+m;
  
  double obj_val,obj_val1,obj_val2,obj_val3;
  arma::mat I, I1, I2, I3;
  double det,det1,det2,det3;
  Rcpp::NumericMatrix para_inf,para_inf1,para_inf2,para_inf3;
  
  // Bayesian estimation under SEL
  TK_Obj rb(X,R,N,-1,L,U,T,GM,type,"",1);
  Roptim<TK_Obj> opt("Nelder-Mead");
  opt.control.trace = 0;
  opt.set_hessian(true);
  opt.minimize(rb, para);
  obj_val=-opt.value();
  I=arma::inv(as<arma::mat>(wrap(opt.hessian())));
  para_inf=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I));
  det =arma::det(I);
  
  TK_Obj rb_nu_sel(X,R,N,-1,L,U,T,GM,type,"GEL",1);
  Roptim<TK_Obj> opt_nu_sel("Nelder-Mead");  
  opt_nu_sel.control.trace = 0;
  opt_nu_sel.set_hessian(true);
  opt_nu_sel.minimize(rb_nu_sel, para);
  obj_val1=-opt_nu_sel.value();
  I1=arma::inv(as<arma::mat>(wrap(opt_nu_sel.hessian())));
  para_inf1=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I1));
  det1=arma::det(I1);
  
  TK_Obj rb_xi_sel(X,R,N,-1,L,U,T,GM,type,"GEL",2);
  Roptim<TK_Obj> opt_xi_sel("Nelder-Mead");
  opt_xi_sel.control.trace = 0;
  opt_xi_sel.set_hessian(true);
  opt_xi_sel.minimize(rb_xi_sel, para);
  obj_val2=-opt_xi_sel.value();
  I2=arma::inv(as<arma::mat>(wrap(opt_xi_sel.hessian())));
  para_inf2=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I2));
  det2=arma::det(I2);
  
  
  TK_Obj rb_spmk_sel(X,R,N,-1,L,U,T,GM,type,"GEL",3);
  Roptim<TK_Obj> opt_spmk_sel("Nelder-Mead");
  opt_spmk_sel.control.trace = 0;
  opt_spmk_sel.set_hessian(true);
  opt_spmk_sel.minimize(rb_spmk_sel, para);
  obj_val3=-opt_spmk_sel.value();
  I3=arma::inv(as<arma::mat>(wrap(opt_spmk_sel.hessian())));
  para_inf3=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I3));
  det3=arma::det(I3);
  
  est_sel={std::sqrt(det1/det)*std::exp(N*(obj_val1-obj_val)),std::sqrt(det2/det)*std::exp(N*(obj_val2-obj_val)),std::sqrt(det3/det)*std::exp(N*(obj_val3-obj_val))};
  bias_sel={est_sel(0)-True_para(0),est_sel(1)-True_para(1),est_sel(2)-True_para(2)};
  // Bayesian estimation under GEL

  TK_Obj rb_nu_gel(X,R,N,q,L,U,T,GM,type,"GEL",1);
  Roptim<TK_Obj> opt_nu_gel("Nelder-Mead");
  opt_nu_gel.control.trace = 0;
  opt_nu_gel.set_hessian(true);
  opt_nu_gel.minimize(rb_nu_gel,para);
  obj_val1=-opt_nu_gel.value();
  I1=arma::inv(as<arma::mat>(wrap(opt_nu_gel.hessian())));
  para_inf1=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I1));
  det1=arma::det(I1);
  
  TK_Obj rb_xi_gel(X,R,N,q,L,U,T,GM,type,"GEL",2);
  Roptim<TK_Obj> opt_xi_gel("Nelder-Mead");
  opt_xi_gel.control.trace = 0;
  opt_xi_gel.set_hessian(true);
  opt_xi_gel.minimize(rb_xi_gel, para);
  obj_val2=-opt_xi_gel.value();
  I2=arma::inv(as<arma::mat>(wrap(opt_xi_gel.hessian())));
  para_inf2=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I2));
  det2=arma::det(I2);
  
  TK_Obj rb_spmk_gel(X,R,N,q,L,U,T,GM,type,"GEL",3);
  Roptim<TK_Obj> opt_spmk_gel("Nelder-Mead");
  opt_spmk_gel.control.trace = 0;
  opt_spmk_gel.set_hessian(true);
  opt_spmk_gel.minimize(rb_spmk_gel, para);
  obj_val3=-opt_spmk_gel.value();
  I3=arma::inv(as<arma::mat>(wrap(opt_spmk_gel.hessian())));
  para_inf3=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I3));
  det3=arma::det(I3);
  est_gel={pow(std::sqrt(det1/det)*std::exp(N*(obj_val1-obj_val)),-1.0/q),pow(std::sqrt(det2/det)*std::exp(N*(obj_val2-obj_val)),-1.0/q),pow(std::sqrt(det3/det)*std::exp(N*(obj_val3-obj_val)),-1.0/q)};
  bias_gel={est_gel(0)-True_para(0),est_gel(1)-True_para(1),est_gel(2)-True_para(2)};
  
  // Bayesian estimation under linex
  TK_Obj rb_nu_lin(X,R,N,c,L,U,T,GM,type,"LINEX",1);
  Roptim<TK_Obj> opt_nu_lin("Nelder-Mead");
  opt_nu_lin.control.trace = 0;
  opt_nu_lin.set_hessian(true);
  opt_nu_lin.minimize(rb_nu_lin, para);
  obj_val1=-opt_nu_lin.value();
  I1=arma::inv(as<arma::mat>(wrap(opt_nu_lin.hessian())));
  para_inf1=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I1));
  det1=arma::det(I1);
  
  TK_Obj rb_xi_lin(X,R,N,c,L,U,T,GM,type,"LINEX",2);
  Roptim<TK_Obj> opt_xi_lin("Nelder-Mead");
  opt_xi_lin.control.trace = 0;
  opt_xi_lin.set_hessian(true);
  opt_xi_lin.minimize(rb_xi_lin, para);
  obj_val2=-opt_xi_lin.value();
  I2=arma::inv(as<arma::mat>(wrap(opt_xi_lin.hessian())));
  para_inf2=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I2));
  det2=arma::det(I2);
  
  TK_Obj rb_spmk_lin(X,R,N,c,L,U,T,GM,type,"LINEX",3);
  Roptim<TK_Obj> opt_spmk_lin("Nelder-Mead");
  opt_spmk_lin.control.trace = 0;
  opt_spmk_lin.set_hessian(true);
  opt_spmk_lin.minimize(rb_spmk_lin, para);
  obj_val3=-opt_spmk_lin.value();
  I3=arma::inv(as<arma::mat>(wrap(opt_spmk_lin.hessian())));
  para_inf3=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I3));
  det3=arma::det(I3);
  
  est_lin={(-1.0/c)*std::log(std::sqrt(det1/det)*std::exp(N*(obj_val1-obj_val))),(-1.0/c)*log(std::sqrt(det2/det)*std::exp(N*(obj_val2-obj_val))),(-1.0/c)*std::log(std::sqrt(det3/det)*std::exp(N*(obj_val3-obj_val)))};
  bias_lin={est_lin(0)-True_para(0),est_lin(1)-True_para(1),est_lin(2)-True_para(2)};
  return Rcpp::List::create(
    Rcpp::Named("est_sel") = est_sel,
    Rcpp::Named("est_gel") = est_gel,
    Rcpp::Named("est_lin") = est_lin,
    Rcpp::Named("bias_sel") = bias_sel,
    Rcpp::Named("bias_gel") = bias_gel,
    Rcpp::Named("bias_lin") = bias_lin);
};

