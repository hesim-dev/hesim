// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
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


// [[Rcpp::export]]
std::vector<arma::mat> list_to_vecmats(Rcpp::List L){
  int n_mats = L.size();
  std::vector<arma::mat> x(n_mats);
  for (int i = 0; i < n_mats; ++i){
    x[i] = Rcpp::as<arma::mat>(L[i]);
  }
  return(x);
}

// [[Rcpp::export]]
std::vector<double> list_loop(Rcpp::List L, arma::rowvec beta, int n){
//   std::vector<double> y;
//   for (int i = 0; i < n; i++){
//    //arma::mat x = Rcpp::as<arma::mat>(L[0]);
// arma::mat x(
//     REAL(VECTOR_ELT(L, 0)), Rf_nrows(VECTOR_ELT(L, 0)), Rf_ncols(VECTOR_ELT(L, 0)), false, true
// );
//     double xb = dot(x.row(0), beta);
  std::vector<arma::mat> x = list_to_vecmats(L);
  std::vector<double> y;
  for (int i = 0; i < n; i++){
    double xb = dot(x[0].row(0), beta);
    y.push_back(xb);
  }
  return(y);
}

// [[Rcpp::export]]
std::vector<double> cube_loop(arma::cube x, arma::rowvec beta, int n){
  std::vector<double> y;
  for (int i=0; i < n; i++){;
    double xb = dot(x.slice(0).row(0), beta);
    y.push_back(xb);
  }
  return(y);
}

/*** R
n <- 100000000

# list
l <- list(matrix(c(1, 2, 3, 4), nrow = 2))
beta <- c(1, 1)
ptm <- proc.time()
tmp <- hesim:::list_loop(l, beta, n)
proc.time() - ptm

# matrix
x <- array(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), dim = c(4,4,1))
beta <- c(1, 1, 1, 1)
ptm <- proc.time()
tmp <- hesim:::cube_loop(x, beta, n)
proc.time() - ptm
*/