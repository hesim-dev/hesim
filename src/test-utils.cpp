// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/utils.h>

/** @defgroup test Tests
 *  Tests of C++ code. These functions are exported to @c R using the
 *  @c R package @c Rcpp and the @c R package @c testthat.
 */

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