// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;

// Convert vector to matrix by row
// [[Rcpp::export]]
arma::mat matrix_byrow(arma::rowvec v, int nrow, int ncol){
  int l = v.n_elem;
  arma::mat m1(0, l);
  m1.insert_rows(0, v);
  m1.reshape(ncol, nrow);
  m1 = arma::trans(m1);
  return(m1);
}

// Convert vector to matrix by column
// [[Rcpp::export]]
arma::mat matrix_bycol(arma::rowvec v, int nrow, int ncol){
  int l = v.n_elem;
  arma::mat m1(0, l);
  m1.insert_rows(0, v);
  m1.reshape(nrow, ncol);
  return(m1);
}

// Convert Rcpp list to std::vector<arma::mat> >
vecmats list_to_vecmats(Rcpp::List l){
  vecmats v;
  int nmats = l.size();
  v.reserve(nmats);
  for (int i = 0; i < nmats; ++i){
    v.push_back(as<arma::mat >(l[i]));
  }
  return v;
}

