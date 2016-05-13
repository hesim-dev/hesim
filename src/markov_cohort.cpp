// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat matrixC(arma::vec v, int nrow, int ncol){
  arma::mat m1;
  m1.insert_cols(0, v);
  m1.reshape(nrow, ncol);
  return(m1);
}

// Bayesian Markov Cohort Simulation
// [[Rcpp::export]]
List bayesianMarkovCohortC(arma::rowvec z0, int ncycles, double discount, int nsims,
                           arma::cube P, arma::cube costs, arma::cube qol,
                           arma::vec P_indx, arma::vec cost_indx, arma::vec qol_indx) {
  // Initialize
  double N = sum(z0);
  int n_states = z0.n_elem;
  int P_len = P.n_slices;
  arma::cube Z(ncycles + 1, n_states, nsims);
  arma::cube Pmat(n_states, n_states, P_len);
  arma::mat cycle_costs(ncycles + 1, nsims);
  arma::mat cycle_qalys(ncycles + 1, nsims);

  // Discount rate
  arma::vec delta(ncycles + 1);
  for (int t = 0; t <= ncycles; ++t){
    delta(t) = pow(1/(1 + discount), t);
  }

  // SIMULATION
  for (int s = 0; s < nsims; ++s){
    // Cycle 0
    Z.slice(s).row(0) = z0;
    cycle_costs(0, s) = 0;
    cycle_qalys(0, s) = 0;

    // Convert vectors to matrices for each unique matrix for simulation s
    for (int j = 0; j < P_len; ++j){
      Pmat.slice(j) = matrixC(P.slice(j).col(s), n_states, n_states);
    }

    // Iterate over cycles 1 - T
    for (int t = 1; t <= ncycles; ++t){
      Z.slice(s).row(t) = Z.slice(s).row(t-1) * Pmat.slice(P_indx(t -1) - 1);
      cycle_costs(t, s) = delta(t) * dot(Z.slice(s).row(t),
                  costs.slice(cost_indx(t - 1) - 1).col(s))/N;
      cycle_qalys(t, s) = delta(t) * dot(Z.slice(s).row(t),
                  qol.slice(qol_indx(t - 1) - 1).col(s))/N;
    }
  }
  return List::create(Named("costs") = cycle_costs,
                      Named("qalys") = cycle_qalys,
                      Named("state") = Z);
}

