// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;

// Predicted probabilities for single individual from multinomial logistic regression
//' @export
// [[Rcpp::export]]
arma::rowvec mlogit_prob(arma::mat beta, arma::rowvec x, int ncat) {
  int zero = 0;
  arma::rowvec odds(ncat);
  arma::rowvec p;
  for (int j = 0; j < ncat; ++j){
    if (j == 0){
      odds(j) = exp(zero);
    }
    else{
      odds(j) = exp(dot(x, beta.row(j - 1)));
    }
  }
  double odds_sum = sum(odds);
  p = odds/odds_sum;
  return p;
}

// Posterior predictive distribution for multinomial logit fit with MCMCpack
//' @export
// [[Rcpp::export]]
arma::cube predict_MCMCmnlC(arma::mat beta, arma::mat x, int ncat) {
  int n = x.n_rows;
  int nsims = beta.n_rows;
  arma::cube prob(n, ncat , nsims);
  int beta_s_cols = beta.n_cols/(ncat - 1);
  arma::mat beta_s(ncat - 1, beta_s_cols);
  for (int s = 0; s < nsims; ++s){
    for (int i = 0; i < n; ++i){
      arma::mat beta_s = matrix_bycol(beta.row(s), ncat - 1, beta_s_cols);
      prob.slice(s).row(i) =  mlogit_prob(beta_s, x.row(i), ncat);
    }
  }
  return(prob);
}
