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
  arma::mat betamat = matrixC(beta, k, nstates - 1);
  arma::rowvec onevec(1); onevec.ones();
  arma::rowvec odds = join_rows(onevec, exp(x * betamat));
  double sum = accu(odds);
  return odds/sum;
}

// Multinomial logit transition probabilities for one individual
// (one row for each origin state for beta)
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

//Markov Cohort Simulation using Multinomial Logit
//' @export
// [[Rcpp::export]]
arma::mat markov_mlogit(arma::mat x, arma::cube beta, arma::rowvec z0, int ncycles, int maxage) {
  int n_states = beta.n_slices;
  int nsims = beta.slice(0).n_rows;
  int beta_cols = beta.n_cols;
  int N = x.n_rows;
  arma::mat beta_s(n_states, beta_cols);
  arma::mat pmat(n_states, n_states);
  arma::mat Z(ncycles + 1, n_states);
  for (int s = 0; s < nsims; ++s){
    for (int j = 0; j < n_states; ++j){
      beta_s.row(j) = beta.slice(j).row(s);
    }
    for (int i = 0; i < N; ++i){
      Z.row(0) = z0;
      for (int t = 1; t <= ncycles; ++t){
        pmat = mlogit_transprob(x.row(i), beta_s, n_states);
        Z.row(t) = Z.row(t - 1) * pmat;
      }
    }
  }
  return Z;
}

// Markov Cohort Simulation with Time Constant Transition Probability
//' @export
// [[Rcpp::export]]
arma::mat markov_pmatC(arma::cube pmat, arma::vec pmat_index, arma::mat z0,
                      int ncycles, int maxage, int nstates) {
  int nsims = pmat.slice(0).n_rows;
  int ncohorts = z0.n_rows;
  int N = (ncycles + 1) * ncohorts * nsims;
  arma::mat Z(N, nstates);
  arma::mat pmati(nstates, nstates);
  int counter = 0;
  for (int s = 0; s < nsims; ++s){
    for (int i = 0; i < ncohorts; ++i){
      pmati = matrixC(pmat.slice(pmat_index(i)), nstates, nstates);
      for (int t = 0; t <= ncycles; ++t){
        if (t == 0){
          Z.row(counter) = z0.row(i);
        }
        else{
          Z.row(counter) = Z.row(counter - 1) * pmati;
        }
        ++counter;
      }
    }
  }
  return Z;
}

// //Bayesian Markov Cohort Simulation
// //' @export
// // [[Rcpp::export]]
// arma::mat markov_transC(arma::mat z0, int ncycles, arma::mat P) {
//   // Initialize
//   int n_states = z0.n_cols;
//   int nsims = z0.n_rows;
//   int N = (ncycles + 1) * nsims;
//   arma::mat Z(N, n_states);
//
//   // SIMULATION
//   int counter = 0;
//   int counter2 = 0;
//   for (int s = 0; s < nsims; ++s){
//     for (int t = 0; t <= ncycles; ++t){
//       if(t == 0){
//         Z.row(counter) = z0.row(s);
//       }
//       else{
//         arma::mat Pmat = matrixC(P.row(counter2), n_states, n_states);
//        Z.row(counter) = Z.row(counter - 1) * Pmat;
//       }
//       ++counter;
//       if (t > 0){
//         ++counter2;
//       }
//     }
//   }
//   return Z;
// }

