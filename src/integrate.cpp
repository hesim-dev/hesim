#include "integrate.h"

class IntegrateTestFunc{
  public:
    double operator()(double x) const {
      return x * x;
    }
};

// [[Rcpp::export]]
double C_test_trapzfun(std::vector<double> x){
  IntegrateTestFunc fun;
  return trapzfun(fun, x.begin(), x.end());
}

// [[Rcpp::export]]
double C_test_trapz(std::vector<double> x, std::vector<double> y){
  return trapz(x.begin(), x.end(), y.begin());
}

// [[Rcpp::export]]
std::vector<double> C_test_cumtrapzfun(std::vector<double> t){
  IntegrateTestFunc fun;
  return cumtrapzfun(fun, t.begin(), t.end());
}

// [[Rcpp::export]]
double C_test_simpsfun(std::vector<double> t){
  IntegrateTestFunc fun;
  return simpsfun(fun, t.begin(), t.end());
}

// [[Rcpp::export]]
std::vector<double> C_test_cumsimpsfun(std::vector<double> t){
  IntegrateTestFunc fun;
  return cumsimpsfun(fun, t.begin(), t.end());
}