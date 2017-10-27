// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;

// Calculate index of maximum of vector
// [[Rcpp::export]]
int vecmax_index(std::vector<double> x) {
  std::vector<double>::iterator it = std::max_element(x.begin(), x.end());
  return it - x.begin();
}

// Calculate maximum of vector
// [[Rcpp::export]]
double vecmax(std::vector<double> x) {
  std::vector<double>::iterator it = std::max_element(x.begin(), x.end());
  return *it;
}

// row maximum of matrix
// [[Rcpp::export]]
arma::colvec rowmaxC(arma::mat x) {
  return arma::max(x, 1);
}

// test mean
// [[Rcpp::export]]
double stdmean(std::vector<double> v) {
  double mean = accumulate(v.begin(), v.end(), 0.0)/v.size();
  return mean;
}

// row maximum of matrix
// [[Rcpp::export]]
arma::ucolvec rowmax_indC(arma::mat x) {
  return arma::index_max(x,1);;
}

// Incremental effect in piecewise comparisons
// [[Rcpp::export]]
std::vector<double> incr_effectC(std::vector<double> x, std::vector<double> y,
                                int nsims, int nstrategies, int ngrps){
  int N = x.size();
  std::vector<double> incr_vec;
  incr_vec.reserve(N);

  // estimate incremental cost and qalys
  double incr = 0;
  int counter = 0;
  for (int g = 0; g < ngrps; ++g){
    for (int j = 0; j < nstrategies; ++j){
      for (int s = 0; s < nsims; ++s){
        incr = x[counter] - y[g * nsims + s];
        incr_vec.push_back(incr);
        ++counter;
      }
    }
  }
  return incr_vec;
}


// Probability CE in pairwise comparisons
// [[Rcpp::export]]
std::vector<double> ceacC(std::vector<double> k, std::vector<double> ie,
                            std::vector<double> ic, int nsims, int nstrategies, int ngrps) {

  // Initialize
  int n_k = k.size();
  int sumpos = 0;
  std::vector<double> prob;
  prob.reserve(n_k * nstrategies * ngrps);
  std::vector<double> k_vec;
  k_vec.reserve(n_k * nstrategies);
  double nb = 0;

  // loop over willingess to pay
  for (int j = 0; j < n_k; ++j){
    int counter = 0;

    // loop over groups
    for (int g = 0; g < ngrps; ++ g){

      // loop over treatment strategies
      for (int i = 0; i < nstrategies; ++i){
        sumpos = 0;

        // loop over simulations
        for (int s = 0; s < nsims; ++s){
           nb = k[j] * ie[counter] - ic[counter];
          if (nb > 0.0){
            ++sumpos;
          }
          ++ counter;
        }
        prob.push_back((double)sumpos/nsims);
        k_vec.push_back(k[j]);
      }
    }
  }
  return prob;
}

// Probability most cost-effective therapy
// [[Rcpp::export]]
std::vector<double> mceC(std::vector<double> k, std::vector<double> e,
                 std::vector<double> c, int nsims, int nstrategies, int ngrps) {

  int n_k = k.size();
  int N = n_k * nstrategies * ngrps;
  std::vector<double> prob(N, 0.0);
  int jg = 0;

  // loop over willingess to pay k
  for (int j = 0; j < n_k; ++j){
    int sg = 0;

    // loop over groups
    for (int g = 0; g < ngrps; ++g){

      // loop over simulations
      for (int s = 0; s < nsims; ++s){
        std::vector<double> nb;
        nb.reserve(nstrategies);

        // loop over treatment strategies
        for (int i = 0; i < nstrategies; ++i){
          nb.push_back(k[j] * e[sg * nstrategies + i] - c[sg * nstrategies + i]);
        }
        int nb_ind = vecmax_index(nb);
        prob[jg * nstrategies + nb_ind] = prob[jg * nstrategies + nb_ind] + 1;
        ++sg;
      }
      ++jg;
    }
  }

  // Convert sum to proportion
  for (int i = 0; i < N; ++i){
    prob[i] = prob[i]/nsims;
  }

  return prob;
}

// Calculate Vstar - expected benefits if choosing optimal arm at each simulation
// [[Rcpp::export]]
std::vector<double> VstarC(std::vector<double> k,
                          std::vector<double> e, std::vector<double> c,
                          int nsims, int nstrategies, int ngrps) {

  int n_k = k.size();
  std::vector<double> Vstar;
  Vstar.reserve(n_k * ngrps);

  // loop over willingess to pay k
  for (int j = 0; j < n_k; ++j){
    int sg = 0;

    // loop over groups
    for (int g = 0; g < ngrps; ++g){
      std::vector<double> Ustar;
      Ustar.reserve(nsims);

      // loop over simulations
      for (int s = 0; s < nsims; ++s){
        std::vector<double> inb;
        inb.reserve(nstrategies);

        // loop over treatment strategies
        for (int i = 0; i < nstrategies; ++i){
          inb.push_back(k[j] * e[sg * nstrategies + i] - c[sg * nstrategies + i]);
        }
        Ustar.push_back(vecmax(inb));
        ++sg;
      }
      double Ustar_mean = accumulate(Ustar.begin(), Ustar.end(), 0.0)/Ustar.size();
      Vstar.push_back(Ustar_mean);
    }
  }

  return Vstar;
}

