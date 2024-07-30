#ifndef LIBFUN_H
#define LIBFUN_H

#include "gamma_fun.h"
#include "nakagami.h"

// The likelihood, ps and posteriors functions and their logarithms 
Nakagami nakagami;

// [[Rcpp::export]]
arma::vec arma_sort(arma::vec x) {
  NumericVector xx=Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(x));
  NumericVector y = clone(xx);
  std::sort(y.begin(), y.end());
  return as<arma::vec>(wrap(y));
}

// [[Rcpp::export]]
double MSE(double eta,Rcpp::NumericVector  eta_hat, std::string type, double q){
  int n=eta_hat.size();
  Rcpp::NumericVector mse(n);
  for(int i=0;i<n;i++){
    if(type=="SEL") mse(i)=pow(eta_hat(i)-eta,2);
    if(type=="GEL") mse(i)= pow(eta_hat(i)/eta,q)-q*std::log(eta_hat(i)/eta)-1;
    if(type=="LINEX") mse(i)=std::exp(q*(eta_hat(i)-eta))-q*(eta_hat(i)-eta)-1;
  };
  return sum(mse)/n;
}

double prior(Rcpp::NumericVector para){
  return std::sqrt(para[0]*gammafun.my_trigamma(para[0])-1)/para[1];
} 

// These classes are used to compute derivative of incomplete gamma
class P1: public Func {
private:
  double nu;
  double xi;
public:
  P1(double a_, double b_) : nu(a_), xi(b_) {}
  
  double operator()(const double& x) const {
    return std::log(x)*std::pow(x,nu-1)*exp(-x);
  };
};

class P2: public Func {
private:
  double nu;
  double xi;
public:
  P2(double a_, double b_) : nu(a_), xi(b_) {}
  
  double operator()(const double& x) const {
    return std::log(x)*std::log(x)*std::pow(x,nu-1)*exp(-x);
  };
};

// The definition of S'_pmk function
// [[Rcpp::export]]
double spmk_fun(Rcpp::NumericVector para, double l, double u, double t, double gm){
  return(R::qnorm(1.0-(1.0-(nakagami.cdf(u,para)-nakagami.cdf(l,para)))/2,0.0,1.0,1,0)/(3.0*std::sqrt(1.0+(2/(pow(para[1],2)*pow(gm,2)))*(std::exp(gm*(para[0]-t))-gm*(para[0]-t)-1))));
}


  
// The mean of ND
// [[Rcpp::export]]
double spmk_mu(double x, double y){
  return (gammafun.my_gamma(x+0.5)/gammafun.my_gamma(x))*std::sqrt(y/x);
}

// The standard deviation of ND
// [[Rcpp::export]]
double spmk_sig(double x, double y){
  return(std::sqrt(y-std::pow((gammafun.my_gamma(x+0.5)/gammafun.my_gamma(x))*std::sqrt(y/x),2)));
};

// [[Rcpp::export]]
Rcpp::NumericVector spmk_grad(Rcpp::NumericVector para, double l, double u, double t, double gm){
  
  double nu=para[0];
  double xi=para[1];
  
  double err_est1;
  int err_code1;
  double err_est2;
  int err_code2;
  
  P1 f1(nu, xi);
  
  double A=std::exp(gm*(nu-t))-gm*(nu-t)-1.0;
  double B=1.0+2.0*A/(pow(gm,2.0)*pow(xi,2.0));
  double p=1.0-(nakagami.cdf(u,{xi,nu})-nakagami.cdf(l,{xi,nu}));
  double Q=R::qnorm(1.0-p/2.0,0.0,1.0,1.0,0);
  double q=1.0/R::dnorm(R::qnorm(1.0-p/2.0,0.0,1.0,1.0,0.0),0.0,1.0,0.0);
  double I_l = integrate(f1, nu*l*l/xi, 1000, err_est1, err_code1);
  double I_u = integrate(f1, nu*u*u/xi, 1000, err_est2, err_code2);
  double psi_l=std::log(gammafun.my_incgam(nu*pow(l,2)/xi,nu));
  double psi_u=std::log(gammafun.my_incgam(nu*pow(u,2)/xi,nu));
  double psi_nu_l=(-pow(nu,nu-1)*std::exp(-nu*l*l/xi)*pow(l,2*nu)*pow(xi,-nu)+I_l)/gammafun.my_incgam(nu*l*l/xi,nu);
  double psi_nu_u=(-pow(nu,nu-1)*std::exp(-nu*u*u/xi)*pow(u,2*nu)*pow(xi,-nu)+I_u)/gammafun.my_incgam(nu*u*u/xi,nu);
  double psi_xi_l=pow(nu,nu)*std::exp(-nu*l*l/xi)*pow(l,2*nu)*pow(xi,-nu-1)/gammafun.my_incgam(nu*l*l/xi,nu);
  double psi_xi_u=pow(nu,nu)*std::exp(-nu*u*u/xi)*pow(u,2*nu)*pow(xi,-nu-1)/gammafun.my_incgam(nu*u*u/xi,nu);
  
  double g1= q/(6.0*std::sqrt(B));
  double g2=(gammafun.my_incgam(nu*l*l/xi,nu)*(psi_nu_l-gammafun.my_digamma(nu))-gammafun.my_incgam(nu*u*u/xi,nu)*(psi_nu_u-gammafun.my_digamma(nu)))/gammafun.my_gamma(nu);
  double g3=(gammafun.my_incgam(nu*l*l/xi,nu)*psi_xi_l-gammafun.my_incgam(nu*u*u/xi,nu)*psi_xi_u)/gammafun.my_gamma(nu);
  double GG1=-(Q/3.0)*(std::exp(gm*(nu-t))-1)/(gm*pow(xi,2)*pow(B,3.0/2.0));
  double GG2= (2.0/3.0)*Q*(std::exp(gm*(nu-t))-gm*(nu-t)-1.0)/(pow(gm,2.0)*pow(xi,3.0)*pow(B,3.0/2.0));  
  double G1=g1*g2+GG1;
  double G2=g1*g3+GG2;
  return {G1,G2};
};

// These class will be used in computing the expectations 
// in EM and informtion using the missing information principle
// expectation E(X)
class E0: public Func {
private:
  double nu;
  double xi;
public:
  E0(double a_, double b_) : nu(a_), xi(b_) {}
  
  double operator()(const double& X) const {
    return std::pow(X,2.0*nu-1.0)*std::exp(-nu*std::pow(X,2.0)/xi);
  };
};

// expectation E(log(X))
class E1: public Func {
private:
  double nu;
  double xi;
public:
  E1(double a_, double b_) : nu(a_), xi(b_) {}
  
  double operator()(const double& x) const {
    return std::log(x)*pow(x,2.0*nu-1.0)*std::exp(-nu*std::pow(x,2.0)/xi);
  };
};

// expectation E(X^2)
class E2: public Func {
private:
  double nu;
  double xi;
public:
  E2(double a_, double b_) : nu(a_), xi(b_) {}
  
  double operator()(const double& x) const {
    return std::pow(x,2)*pow(x,2.0*nu-1.0)*exp(-nu*pow(x,2.0)/xi);
  };
};

// expectation E(log(X)^2)
class E3: public Func {
private:
  double nu;
  double xi;
public:
  E3(double a_, double b_) : nu(a_), xi(b_) {}
  
  double operator()(const double& x) const {
    return std::log(x)*std::log(x)*pow(x,2.0*nu-1.0)*exp(-nu*pow(x,2.0)/xi);
  };
};

// expectation E(X^4)
class E4: public Func {
private:
  double nu;
  double xi;
public:
  E4(double a_, double b_) : nu(a_), xi(b_) {}
  
  double operator()(const double& x) const {
    return std::pow(x,4)*pow(x,2.0*nu-1.0)*exp(-nu*pow(x,2.0)/xi);
  };
};

// expectation E(X^2*log(X))
class E5: public Func {
private:
  double nu;
  double xi;
public:
  E5(double a_, double b_) : nu(a_), xi(b_) {}
  
  double operator()(const double& x) const {
    return std::log(x)*std::pow(x,2)*pow(x,2.0*nu-1.0)*exp(-nu*pow(x,2.0)/xi);
  };
};

// likelihood function
// [[Rcpp::export]]
double like(Rcpp::NumericVector para, 
            Rcpp::NumericVector X, 
            Rcpp::NumericVector R){
  
  // Compute objective value
  double lk=1;
  int m=X.size();
  for(int i=0;i<m;i++){
    lk *= nakagami.pdf(X(i),para)*pow(nakagami.sur(X(i),para),R(i));
  }
  return lk;
}  

// loglikelihood function
// [[Rcpp::export]]
double loglike(Rcpp::NumericVector para, 
            Rcpp::NumericVector X, 
            Rcpp::NumericVector R){
  
  // Compute objective value
  double lk=1;
  int m=X.size();
  for(int i=0;i<m;i++){
    lk *= nakagami.pdf(X(i),para)*pow(nakagami.sur(X(i),para),R(i));
  }
  return std::log(lk);
}  

// [[Rcpp::export]]
double logpostlk(Rcpp::NumericVector para, 
            Rcpp::NumericVector X, 
            Rcpp::NumericVector R){
  
  // Compute objective value
  double lk=1;
  int m=X.size();
  for(int i=0;i<m;i++){
    lk *= nakagami.pdf(X(i),para)*pow(nakagami.sur(X(i),para),R(i));
  }
  return std::log(lk)+std::log(prior(para));
}  

// [[Rcpp::export]]
double mps(Rcpp::NumericVector para, 
           Rcpp::NumericVector X, 
           Rcpp::NumericVector R){
  
  // Compute objective value
  int m=X.size();
  double mp=nakagami.sur(X(m-1),para)*nakagami.cdf(X(0),para)*pow(nakagami.sur(X(0),para),R(0));
  for(int i=1;i<m;i++){
    mp *= (nakagami.cdf(X(i),para)-nakagami.cdf(X(i-1),para))*pow(nakagami.sur(X(i),para),R(i));
  }
  return mp;
}

// [[Rcpp::export]]
double logmps(Rcpp::NumericVector para, 
           Rcpp::NumericVector X, 
           Rcpp::NumericVector R){
  
  // Compute objective value
  int m=X.size();
  double mp=nakagami.sur(X(m-1),para)*nakagami.cdf(X(0),para)*pow(nakagami.sur(X(0),para),R(0));
  for(int i=1;i<m;i++){
    mp *= (nakagami.cdf(X(i),para)-nakagami.cdf(X(i-1),para))*pow(nakagami.sur(X(i),para),R(i));
  }
  return std::log(mp);
}

// [[Rcpp::export]]
double logpostps(Rcpp::NumericVector para, 
           Rcpp::NumericVector X, 
           Rcpp::NumericVector R){
  
  // Compute objective value
  int m=X.size();
  double mp=nakagami.sur(X(m-1),para)*nakagami.cdf(X(0),para)*pow(nakagami.sur(X(0),para),R(0));
  for(int i=1;i<m;i++){
    mp *= (nakagami.cdf(X(i),para)-nakagami.cdf(X(i-1),para))*pow(nakagami.sur(X(i),para),R(i));
  }
  return std::log(mp)+std::log(prior(para));
}

class LogObjective : public Functor {
  private: arma::vec X;
    arma::vec R;
    std::string Type;
public:
  LogObjective(arma::vec xx_, arma::vec rr_, std::string type_) : X(xx_), R(rr_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LK")
      return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R))));
    if(Type=="PS")
      return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R))));
  }
};


class Posterior : public Functor {
  private: arma::vec X;
    arma::vec R;
    std::string Type;
public:
  Posterior(arma::vec xx_, arma::vec rr_, std::string type_) : X(xx_), R(rr_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LK")
      return like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)));
    if(Type=="PS")
      return mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)));;
  }
};

class LogPosterior : public Functor {
  private: arma::vec X;
    arma::vec R;
    std::string Type;
public:
  LogPosterior(arma::vec xx_, arma::vec rr_, std::string type_) : X(xx_), R(rr_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LK")
      return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                  Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                  Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))));
    if(Type=="PS")
      return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                 Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                 Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))));
  }
};

// [[Rcpp::export]]
NumericVector my_as_mcmc(NumericVector x){
  Rcpp::Function as_mcmc_rcpp("as.mcmc");
  Rcpp::NumericVector res=as_mcmc_rcpp(x);
  return res;
}

// [[Rcpp::export]]
NumericVector my_HPD(NumericVector x){
  Rcpp::Function HPD_rcpp("HPDinterval");
  Rcpp::NumericVector res=HPD_rcpp(x);
  return res;
}


class TK_Obj : public Functor {
  private: arma::vec X;
    arma::vec R;
    int N;
    double Q;
    double L;
    double U;
    double T;
    double GM;
    std::string Type;
    std::string Loss;
    int NoPara;
public:
  TK_Obj(arma::vec xx_, arma::vec rr_, int nn_,double qq_, double ll_, double uu_,double tt_,double gm_, std::string type_,std::string loss_, int nopara_) : X(xx_), R(rr_), N(nn_), Q(qq_), L(ll_),U(uu_),T(tt_),GM(gm_),Type(type_), Loss(loss_), NoPara(nopara_){}

  double operator()(const arma::vec &y) override {
    if((Type=="LK") && (Loss==""))
      return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))/N;
    if((Type=="PS") && (Loss==""))
      return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))/N;
    if((Type=="LK")&&(Loss=="GEL")&&(NoPara==1))
      return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*pow(y[0],-Q))/N;
    if((Type=="PS")&&(Loss=="GEL")&&(NoPara==1))
      return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*pow(y[0],-Q))/N;
    if((Type=="LK")&&(Loss=="LINEX")&&(NoPara==1))
      return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*exp(-Q*y[0]))/N;
    if((Type=="PS")&&(Loss=="LINEX")&&(NoPara==1))
      return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*exp(-Q*y[0]))/N;
    if((Type=="LK")&&(Loss=="GEL")&&(NoPara==2))
      return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*pow(y[1],-Q))/N;
    if((Type=="PS")&&(Loss=="GEL")&&(NoPara==2))
      return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*pow(y[1],-Q))/N;
    if((Type=="LK")&&(Loss=="LINEX")&&(NoPara==2))
      return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                            Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*exp(-Q*y[1]))/N;
    if((Type=="PS")&&(Loss=="LINEX")&&(NoPara==2))
      return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                           Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*exp(-Q*y[1]))/N;
      if((Type=="LK")&&(Loss=="GEL")&&(NoPara==3))
        return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                              Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                              Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*pow(spmk_fun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),L,U,T,GM),-Q))/N;
      if((Type=="PS")&&(Loss=="GEL")&&(NoPara==3))
        return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                             Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                             Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*pow(spmk_fun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),L,U,T,GM),-Q))/N;
      if((Type=="LK")&&(Loss=="LINEX")&&(NoPara==3))
        return -std::log(like(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                              Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                              Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*std::exp(-Q*spmk_fun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),L,U,T,GM)))/N;
      if((Type=="PS")&&(Loss=="LINEX")&&(NoPara==3))
        return -std::log(mps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                             Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                             Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)))*prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))*std::exp(-Q*spmk_fun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),L,U,T,GM)))/N;
    };
};

double cloglike(Rcpp::NumericVector para,
                Rcpp::NumericVector X, 
                Rcpp::NumericVector R,
                Rcpp::NumericVector Z1,
                Rcpp::NumericVector Z2){
  double nu=para[0];
  double xi=para[1];
  int m=R.size();
  int n=m+sum(R);
  Rcpp::NumericVector ll(m);
  for(int i=0;i<m;i++)
     ll=(2.0*nu-1.0)*std::log(X[i])-(nu/xi)*pow(X[i],2)+(2.0*nu-1.0)*R[i]*Z1[i]-(nu/xi)*R[i]*Z2[i];
  double cllk=-n*log(gammafun.my_gamma(nu))+n*nu*(log(nu/xi))+sum(ll);
  return -cllk;
} 

class CLogLike : public Functor {
  private: arma::vec X;
    arma::vec R;
    arma::vec Z1;
    arma::vec Z2;
    int n;
public:
  CLogLike(arma::vec xx_, arma::vec rr_,arma::vec z1_,arma::vec z2_,int nn_) : X(xx_), R(rr_), Z1(z1_), Z2(z2_), n(nn_){}
  
  double operator()(const arma::vec &y) override {
    return cloglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)), 
                    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(Z1)),
                    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(Z2)));
  };
};

// Information matrix using missing information principle
// [[Rcpp::export]]
NumericMatrix inform(Rcpp::NumericVector para, 
                     Rcpp::NumericVector X, Rcpp::NumericVector R, 
                     double lower, double upper){
  double nu=para(0);
  double xi=para(1);
  int m=R.size();
  int n=m+sum(R);
  NumericVector W0(m);
  NumericVector W1(m);
  NumericVector W2(m);
  NumericVector W3(m);
  NumericVector W4(m);
  NumericVector W5(m);
  NumericVector Z1(m);
  NumericVector Z2(m);
  NumericVector Z3(m);
  NumericVector Z4(m);
  NumericVector Z5(m);
  NumericVector I1(m);
  
  NumericMatrix hess(2,2);
  
  double err_est0, err_est1, err_est2, err_est3, err_est4, err_est5, err_est6;
  
  int err_code0, err_code1, err_code2, err_code3, err_code4, err_code5, err_code6;
  
  E0 f0(nu, xi);
  E1 f1(nu, xi);
  E2 f2(nu, xi);
  E3 f3(nu, xi);
  E4 f4(nu, xi);
  E5 f5(nu, xi);
  P1 f6(nu, xi);
  
  for(int i=0;i<m;i++){
    I1(i) = integrate(f6, lower, nu*pow(X(i),2)/xi, err_est6, err_code6);
    W0(i) = integrate(f0, X(i), upper, err_est0, err_code0);
    W1(i) = integrate(f1, X(i), upper, err_est1, err_code1);
    W2(i) = integrate(f2, X(i), upper, err_est2, err_code2);
    W3(i) = integrate(f3, X(i), upper, err_est3, err_code3);
    W4(i) = integrate(f4, X(i), upper, err_est4, err_code4);
    W5(i) = integrate(f5, X(i), upper, err_est5, err_code5);
    Z1(i)=W1(i)/W0(i);
    Z2(i)=W2(i)/W0(i);
    Z3(i)=W3(i)/W0(i);
    Z4(i)=W4(i)/W0(i);
    Z5(i)=W5(i)/W0(i);
  };
  
  double a11=n*R::trigamma(nu)-n/nu;
  double a12=0;
  double a22=n*nu/pow(xi,2);
  
  NumericVector psi(m), psi_nu(m), psi_xi(m), h1(m), h2(m);
  for(int i=0;i<m;i++){
    psi(i)=std::log(gammafun.my_incgam(nu*pow(X(i),2)/xi,nu));
    psi_nu(i)=(gammafun.my_digamma(nu)-pow(nu,nu-1)*std::exp(-nu*pow(X(i),2)/xi)*pow(X(i),2*nu)*pow(xi,-nu)-I1(i))/gammafun.my_incgam(nu*pow(X(i),2)/xi,nu);
    psi_xi(i)=pow(nu,nu)*std::exp(-nu*pow(X(i),2)/xi)*pow(X(i),2*nu)*pow(xi,-nu-1)/gammafun.my_incgam(nu*pow(X(i),2)/xi,nu);
    
    h1(i)=std::log(nu/xi)+1-gammafun.my_digamma(nu)-psi_nu(i);
    h2(i)=-nu/xi-psi_xi(i);
  }
  double b11=sum(R*(pow(h1,2)+4*h1*Z1-2*(h1/xi)*Z2+4*Z3+(1/pow(xi,2))*Z4-(4/xi)*Z5));  
  double b22=sum(R*(pow(h2,2)+2*h2*(nu/pow(xi,2))*Z2+(pow(nu,2)/pow(xi,4))*Z4));  
  double b12=sum(R*(h1*h2+2*h2*Z1+(nu*h1/pow(xi,2)-h2/xi)*Z2-(nu/pow(xi,3))*Z4+(2*nu/pow(xi,2))*Z5));  
  
  hess(0,0)=a11-b11;
  hess(0,1)=a12-b12;
  hess(1,0)=a12-b12;
  hess(1,1)=a22-b22;
  
  //                return hess;
  return hess;
}
// EM algorithm using the solving the normal equation
// [[Rcpp::export]]
double cscore(double y, double xi,
              Rcpp::NumericVector X, 
              Rcpp::NumericVector R,
              Rcpp::NumericVector Z){
  int m=R.size();
  double n=m+sum(R);
  double F=0.0;
  for(int i=0;i<m;i++)
    F+=std::log(X(i))+R(i)*Z(i);
  return gammafun.my_digamma(y)-std::log(y)-2.0*F/n+std::log(xi);
} 

#endif

