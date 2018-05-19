// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/utils.h>

/****************************
* Custom Rcpp::as converters
****************************/
template <typename T1, typename T2> 
T1 list_to_vec(Rcpp::List l){
  T1 v;
  int n = l.size();
  v.reserve(n);
  for (int i = 0; i < n; ++i){
    v.push_back(Rcpp::as<T2 >(l[i]));
  }
  return v;
}

// Convert Rcpp list to vecmats - i.e., std::vector<arma::mat> >
vecmats list_to_vecmats(Rcpp::List l){
  vecmats v;
  int nmats = l.size();
  v.reserve(nmats);
  for (int i = 0; i < nmats; ++i){
    v.push_back(Rcpp::as<arma::mat >(l[i]));
  }
  return v;
}

// Convert Rcpp list to vecmats_2d - i.e., std::vector<vecmats>
vecmats_2d list_to_vecmats_2d(Rcpp::List l){
  vecmats_2d v;
  int n = l.size();
  v.reserve(n);
  for (int i = 0; i < n; ++i){
    v.push_back(list_to_vecmats(l[i]));
  }
  return v;
}

namespace Rcpp {
  template <> vecmats as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return list_to_vecmats(l);
  }
  
  template <> vecmats_2d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return list_to_vecmats_2d(l);
  }
  
  template <> vecmats_3d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return list_to_vec<vecmats_3d, vecmats_2d> (l);
  }
  
  template <> vecstrings_2d as(SEXP object) {
    Rcpp::List l = Rcpp::as<Rcpp::List>(object);
    return list_to_vec<vecstrings_2d, vecstrings> (l);
  }
}

/*************************
* Convert vector to matrix
*************************/
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

/**********************************
* Test adding constant to a vector
*********************************/
// [[Rcpp::export]]
std::vector<int> test_C_add_constant(std::vector<int> v, double value){
  return add_constant(v, value);
}


