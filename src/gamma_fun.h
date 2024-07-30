#ifndef GAMMA_FUN_H
#define GAMMA_FUN_H

class GammaFun{
public:
  GammaFun(){};
  double my_gamma(double x){
    Rcpp::Function gam_rcpp("gamma");
    Rcpp::NumericVector res=gam_rcpp(x);
    return res[0];
  }
  double my_incgam(double x, double y){
    Rcpp::Function incgam_rcpp("incgam");
    Rcpp::NumericVector res=incgam_rcpp(x, y);
    return res[0];
  }
  double my_digamma(double x){
    Rcpp::Function digam_rcpp("digamma");
    Rcpp::NumericVector res=digam_rcpp(x);
    return res[0];
  }
  double my_trigamma(double x){
    Rcpp::Function trigam_rcpp("trigamma");
    Rcpp::NumericVector res=trigam_rcpp(x);
    return res[0];
  }
  double my_tetragamma(double x){
    Rcpp::Function tetragam_rcpp("psigamma");
    Rcpp::NumericVector res=tetragam_rcpp(x,2);
    return res[0];
  }
};
#endif

