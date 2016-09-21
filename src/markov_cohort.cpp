// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Convert vector to matrix by row
//' @export
// [[Rcpp::export]]
arma::mat matrixC(arma::rowvec v, int nrow, int ncol){
  int l = v.n_elem;
  arma::mat m1(0, l);
  m1.insert_rows(0, v);
  m1.reshape(ncol, nrow);
  m1 = arma::trans(m1);
  return(m1);
}

// Multinomial logit predicted probabilities
//' @export
// [[Rcpp::export]]
arma::rowvec mlogit_prob(arma::rowvec x, arma::rowvec beta, int nstates) {
  int k = x.n_elem;
  arma::mat betamat = matrixC(beta, k, nstates -1);
  arma::rowvec onevec(1); onevec.ones();
  arma::rowvec odds = join_rows(onevec, exp(x * betamat));
  double sum = accu(odds);
  return odds/sum;
}

// Multinomial logit transition probabilities
//' @export
// [[Rcpp::export]]
arma::mat mlogit_transprob(arma::rowvec x, arma::mat beta, int nstates) {
  arma::mat pmat(nstates, nstates);
  for (int j = 0; j < nstates; ++j){
    pmat.row(j) =  mlogit_prob(x, beta.row(j), nstates);
  }
  return pmat;
}

// Adjust transition probability for mortality probability
//' @export
// [[Rcpp::export]]
arma::mat transprob_addmort(arma::mat p, arma::vec pmort, int nstates) {
  for (int i = 0; i < nstates; ++i){
    p.row(i) = (1 - pmort(i)) * p.row(i);
  }
  p = join_rows(p, pmort);
  return p;
}

//Bayesian Markov Cohort Simulation
//' @export
// [[Rcpp::export]]
arma::mat markov_trans2(arma::mat x, arma::cube beta, int ncycles, int maxage) {
  return x;
}

//Bayesian Markov Cohort Simulation
//' @export
// [[Rcpp::export]]
arma::mat markov_transC(arma::mat z0, int ncycles, arma::mat P) {
  // Initialize
  int n_states = z0.n_cols;
  int nsims = z0.n_rows;
  int N = (ncycles + 1) * nsims;
  arma::mat Z(N, n_states);

  // SIMULATION
  int counter = 0;
  int counter2 = 0;
  for (int s = 0; s < nsims; ++s){
    for (int t = 0; t <= ncycles; ++t){
      if(t == 0){
        Z.row(counter) = z0.row(s);
      }
      else{
        arma::mat Pmat = matrixC(P.row(counter2), n_states, n_states);
       Z.row(counter) = Z.row(counter - 1) * Pmat;
      }
      ++counter;
      if (t > 0){
        ++counter2;
      }
    }
  }
  return Z;
}

