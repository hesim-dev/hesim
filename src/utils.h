# ifndef UTILS_H
# define UTILS_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat matrix_byrow(arma::rowvec v, int nrow, int ncol);
arma::mat matrix_bycol(arma::rowvec v, int nrow, int ncol);

# endif
