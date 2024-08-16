#ifndef NAKAGAMI_H
#define NAKAGAMI_H


#include "gamma_fun.h"

GammaFun gammafun;

// [[Rcpp::export]]
class Nakagami{
public:
  Nakagami(){};
  // The quantile function of Nakagami Dist
  double qnt(double w, double xi,double nu){
    Rcpp::Function qnaka_rcpp("qnaka");
    Rcpp::NumericVector res=qnaka_rcpp(w,xi,nu);
    return res[0];
  }
  // The density function of Nakagami Dist
  double pdf(double w, double xi,double nu){
    Rcpp::Function dnaka_rcpp("dnaka");
    Rcpp::NumericVector res=dnaka_rcpp(w,xi,nu);
    return res[0];
  }
  // The distribution function of Nakagami Dist
  double cdf(double w, double xi,double nu){
    Rcpp::Function pnaka_rcpp("pnaka");
    Rcpp::NumericVector res=pnaka_rcpp(w,xi,nu);
    return res[0];
  }
  // The survival function of Nakagami Dist
  double sur(double w, double xi, double nu){
    Rcpp::Function pnaka_rcpp("pnaka");
    Rcpp::NumericVector res=pnaka_rcpp(w,xi,nu);
    return 1.0-res[0];
  }
  // The mean of Nakagami Dist
  double mu(double x, double y){
    return (gammafun.my_gamma(x+0.5)/gammafun.my_gamma(x))*std::sqrt(y/x);
  }
  // The standard deviation of Nakagami Dist
  double sig(double x, double y){
    return(std::sqrt(y-std::pow((gammafun.my_gamma(x+0.5)/gammafun.my_gamma(x))*std::sqrt(y/x),2)));
  }
};
#endif

