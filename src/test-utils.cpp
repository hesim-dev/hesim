// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/utils.h>

/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<int> C_test_add_constant_int(std::vector<int> v, double value){
  hesim::add_constant(v, value);
  return v;
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<double> C_test_add_constant_double(std::vector<double> v, double value){
  hesim::add_constant(v, value);
  return v;
}