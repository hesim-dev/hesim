// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/utils.h>

/****************************
* Custom Rcpp::as converters
****************************/
namespace Rcpp {
  template <> vecmats as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return hesim::internal::list_to_vec<vecmats, arma::mat>(l);
  }
  
  template <> vecmats_2d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return hesim::internal::list_to_vec<vecmats_2d, vecmats>(l);
  }
  
  template <> vecmats_3d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return hesim::internal::list_to_vec<vecmats_3d, vecmats_2d> (l);
  }
  
  template <> vecstrings_2d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return hesim::internal::list_to_vec<vecstrings_2d, vecstrings> (l);
  }
}

/************************
* Functions exported to R
************************/
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

/******
* Tests
******/
// [[Rcpp::export]]
std::vector<int> C_test_add_constant_int(std::vector<int> v, double value){
  hesim::add_constant(v, value);
  return v;
}

// [[Rcpp::export]]
std::vector<double> C_test_add_constant_double(std::vector<double> v, double value){
  hesim::add_constant(v, value);
  return v;
}


