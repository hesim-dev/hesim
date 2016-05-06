// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List markovCohort(arma::cube P, arma::rowvec z0, int ncycles,
                  arma::mat state_costs, arma::mat state_qol, double discount) {
  // Initialize
  double N = sum(z0);
  int n_states = P.n_cols;
  int len_state_costs = state_costs.n_rows;
  int len_state_qol = state_qol.n_rows;
  int len_P = P.n_slices;
  arma::mat Z1(ncycles, n_states);
  arma::mat Z = join_cols(z0, Z1);
  arma::vec qalys(ncycles + 1);
  arma::vec costs(ncycles + 1);
  double delta = 1;

  // Cycle 0
  qalys(0) = 0;
  costs(0) = 0;

  // Markov Chain
  for (int t = 1; t <= ncycles; ++t){
    delta = pow(1/(1 + discount), t);

    // state vector
    if(len_P == 1){
      Z.row(t) = Z.row(t-1) * P.slice(0);
    }
    else{
      Z.row(t) = Z.row(t-1) * P.slice(t-1);
    }

    // costs
    if (len_state_costs == 1){
      costs(t) = delta * dot(Z.row(t), state_costs.row(0))/N;
    }
    else {
      costs(t) = delta * dot(Z.row(t), state_costs.row(t))/N;
    }

    // qalys
    if (len_state_qol == 1){
      qalys(t) = delta * dot(Z.row(t), state_qol.row(0))/N;
    }
    else {
      qalys(t) = delta * dot(Z.row(t), state_qol.row(t))/N;
    }
  }

  return List::create(Named("state") = Z,
                      Named("costs") = costs,
                      Named("qalys") = qalys,
                      Named("delta") = delta);
}

