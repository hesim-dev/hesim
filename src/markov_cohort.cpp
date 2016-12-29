// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;

// Adjust transition probability for mortality probability
//' @export
// [[Rcpp::export]]
arma::mat transprob_addmort(arma::mat p, arma::rowvec pmort, int nstates,
                            arma::rowvec zerovec, arma::rowvec onescalar) {
  for (int i = 0; i < nstates; ++i){
    p.row(i) = (1 - pmort(i)) * p.row(i);
  }
  p = join_cols(p, zerovec);
  //arma::rowvec zerovec = arma::zeros<arma::rowvec>(pmat_cols);
  //arma::rowvec onescalar = arma::ones<arma::rowvec>(1);
  pmort = join_rows(pmort, onescalar);
  p = join_rows(p, trans(pmort));
  return p;
}

// Markov cohort simulation
//' @export
// [[Rcpp::export]]
arma::mat markov_cohort_transC(arma::mat z0, arma::vec ncycles,
                               arma::cube &pmat, arma::vec pmat_index,
                               bool mortadj, arma::cube mortprob, arma::vec mortprob_index){
  int ncohorts = z0.n_rows;
  int nstates = z0.n_cols;
  int pmat_cols = pmat.slice(0).n_cols;
  arma::mat pmat_i(pmat_cols, pmat_cols);
  int N = sum(ncycles) + ncohorts;
  arma::mat Z(N, nstates);
  arma::rowvec mortprob_it(pmat_cols);
  arma::mat pmat_iadj(nstates, nstates);
  arma::rowvec zerovec = arma::zeros<arma::rowvec>(pmat_cols);
  arma::rowvec onevec = arma::ones<arma::rowvec>(1);
  int counter = 0;
  for (int i = 0; i < ncohorts; ++i){
    pmat_i = pmat.slice(pmat_index(i));
    for (int t = 0; t <= ncycles(i); ++t){
      if (t == 0){
        Z.row(counter) = z0.row(i);
      }
      else{
        if (mortadj){
          mortprob_it = mortprob.slice(mortprob_index(i)).row(t);
          pmat_iadj = transprob_addmort(pmat_i, mortprob_it, pmat_cols,
                                        zerovec, onevec);
          Z.row(counter) = Z.row(counter - 1) * pmat_iadj;
        }
        else{
          Z.row(counter) = Z.row(counter - 1) * pmat_i;
        }
      }
      ++counter;
    }
  }
  return Z;
}

