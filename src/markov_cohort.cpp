// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List markovCohortC(arma::cube P, arma::rowvec z0, int ncycles,
                  arma::mat costs, arma::mat qol, double discount,
                  arma::vec P_indx, arma::vec cost_indx, arma::vec qol_indx) {
  // Initialize
  double N = sum(z0);
  int n_states = P.n_cols;
  arma::mat Z1(ncycles, n_states);
  arma::mat Z = join_cols(z0, Z1);
  arma::vec cycle_qalys(ncycles + 1);
  arma::vec cycle_costs(ncycles + 1);
  double delta = 1;

  // Cycle 0
  cycle_qalys(0) = 0;
  cycle_costs(0) = 0;

  //Cycles 1 - T
  for (int t = 1; t <= ncycles; ++t){
    delta = pow(1/(1 + discount), t);
    Z.row(t) = Z.row(t - 1) * P.slice(P_indx(t - 1) - 1);
    cycle_costs(t) = delta * dot(Z.row(t), costs.row(cost_indx(t - 1) - 1))/N;
    cycle_qalys(t) = delta * dot(Z.row(t), qol.row(qol_indx(t - 1) - 1))/N;
  }

  return List::create(Named("state") = Z,
                      Named("costs") = cycle_costs,
                      Named("qalys") = cycle_qalys);
}

