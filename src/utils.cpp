// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/utils.h>
// #include <gperftools/profiler.h> // Profiling on Linux

/**
 * Calculate the maximum value of each row in a matrix
 * @param x A matrix from the Armadillo library
 * @return A vector of maximum values with length equal to the number
 * of rows in x
 */
// [[Rcpp::export]]
arma::colvec C_rowmax(arma::mat x) {
  return arma::max(x, 1);
}

/**
 * Calculate the indices of maximum values in each row of a matrix
 * @param x A matrix from the Armadillo library
 * @return A vector of the column indices of the maximum values of each row with
 * length equal to the number of rows in x
 */
// [[Rcpp::export]]
arma::ucolvec C_rowmax_index(arma::mat x) {
  return arma::index_max(x,1);
}

// // Profiling on Linux using Google perftools
// // [[Rcpp::export]]
// SEXP start_profiler(std::string str) {
//   ProfilerStart(str.c_str());
//   return R_NilValue;
// }

// // [[Rcpp::export]]
// SEXP stop_profiler() {
//   ProfilerStop();
//   return R_NilValue;
// }
