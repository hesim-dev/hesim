// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Convert vector to matrix by row
//' @export
// [[Rcpp::export]]
arma::mat matrix_byrow(arma::rowvec v, int nrow, int ncol){
  int l = v.n_elem;
  arma::mat m1(0, l);
  m1.insert_rows(0, v);
  m1.reshape(ncol, nrow);
  m1 = arma::trans(m1);
  return(m1);
}

// Convert vector to matrix by row
//' @export
// [[Rcpp::export]]
arma::mat matrix_bycol(arma::rowvec v, int nrow, int ncol){
  int l = v.n_elem;
  arma::mat m1(0, l);
  m1.insert_rows(0, v);
  m1.reshape(nrow, ncol);
  return(m1);
}
