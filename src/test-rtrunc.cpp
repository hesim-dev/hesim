// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/stats/distributions.h>

// [[Rcpp::export]]
double C_test_rtrunc_repeat(double lower, double upper){
  hesim::stats::exponential exp(.75);
  auto f_random = [exp](){ return exp.random(); };
  return hesim::stats::rtrunc_repeat(f_random, lower, upper);
}