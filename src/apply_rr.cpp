// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

arma::rowvec apply_complement(arma::rowvec x, const int complement) {
  double rowsums = 0;
  for (int i = 0; i < x.size(); ++i) {
    if (i != complement) rowsums += x[i];
  }
  if (rowsums > 1) {
    Rcpp::stop("Transition probabilities cannot be negative.");
  } 
  x(complement) = 1 - rowsums;
  return x;
}

void apply_complement(arma::mat &x, const arma::umat complement) {
  for (int i = 0; i < complement.n_rows; ++i) {
    int r = complement(i, 0);
    x.row(r) = apply_complement(x.row(r), complement(i, 1));
  }
}

arma::mat apply_rr(const arma::mat &x, const arma::rowvec rr, const arma::umat index,
              const arma::umat complement){
  arma::mat y = x;
  for (int i = 0; i < index.n_rows; ++ i) { 
    int r = index(i, 0); 
    int s = index(i, 1);
    y(r, s) = x(r, s) * rr(i);
  }
  apply_complement(y, complement);
  return y;
}

/**
 * Apply relative risks to transition probability matrices. See the documentation 
 * for the @c R function @c apply_rr() for more details.
 * @param x A cube where each slice is a square transition probability matrix.
 * @param rr A matrix where each column is a vector of relative risks to apply
 * to each transition matrix in @p x.  
 * @param index A matrix of integers denoting the indices of the matrices in 
 * @p x that @p rr is applied to. The first column denotes a row and the 
 * second column a column.
 * @param complement A matrix of integers in the same format as @p x 
 * denoting elements of the matrices in @p x that are "complements" 
 * (i.e., computed as 1 less the sum of all other elements in that row). There
 * can be at most one complementary column for each row of a matrix @p x. 
 * @return The same as described in the @c R function @c apply_rr().
 */
// [[Rcpp::export]]
arma::cube C_apply_rr(const arma::cube &x, const arma::mat rr, const arma::umat index,
                     const arma::umat complement) {
  int N = rr.n_rows;
  arma::cube y(x.n_rows, x.n_cols, N); 
  for (int i = 0; i < N; ++i) {
    y.slice(i) = apply_rr(x.slice(i % x.n_slices), rr.row(i), index, complement);
  }
  return y;
}

