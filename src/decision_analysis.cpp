// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;

// Calculate index of maximum of vector
//' @export
// [[Rcpp::export]]
int vecmax_index(std::vector<double> x) {
  std::vector<double>::iterator it = std::max_element(x.begin(), x.end());
  return it - x.begin();
}

// Calculate maximum of vector
//' @export
// [[Rcpp::export]]
double vecmax(std::vector<double> x) {
  std::vector<double>::iterator it = std::max_element(x.begin(), x.end());
  return *it;
}

// row maximum of matrix
//' @export
// [[Rcpp::export]]
arma::colvec rowmaxC(arma::mat x) {
  return arma::max(x, 1);
}

// test mean
//' @export
// [[Rcpp::export]]
double stdmean(std::vector<double> v) {
  double mean = accumulate(v.begin(), v.end(), 0.0)/v.size();
  return mean;
}

// row maximum of matrix
//' @export
// [[Rcpp::export]]
arma::ucolvec rowmax_indC(arma::mat x) {
  return arma::index_max(x,1);;
}

// Incremental change in piecewise comparisons
//' @export
// [[Rcpp::export]]
std::vector<double> incr_changeC(std::vector<double> x, std::vector<double> y,
                                int nsims, int narms){
  int N = x.size();
  std::vector<double> incr_vec;
  incr_vec.reserve(N);
  
  // estimate incremental cost and qalys
  double incr = 0;
  int counter = 0;
  for (int j = 0; j < narms; ++j){
    for (int s = 0; s < nsims; ++s){
      incr = x[counter] - y[s];
      incr_vec.push_back(incr);
      ++counter;
    }
  }
  return incr_vec;
}


// Probability CE in pairwise comparisons
//' @export
// [[Rcpp::export]]
std::vector<double> ceacC(std::vector<double> k, std::vector<double> e,
                            std::vector<double> c, std::vector<double> e_comp,
                            std::vector<double> c_comp, int nsims, int narms) {

  // Initialize 
  int n_k = k.size();
  int sumpos = 0;
  std::vector<double> prob;
  prob.reserve(n_k * narms);
  std::vector<double> k_vec;
  k_vec.reserve(n_k * narms);
  double nb = 0;
  
  // Incremental cost and benefit
  std::vector<double> ic = incr_changeC(c, c_comp, nsims, narms);
  std::vector<double> ie = incr_changeC(e, e_comp, nsims, narms);
  
  // loop over willingess to pay
  for (int j = 0; j < n_k; ++j){
    
    // loop over treatment arms
    for (int i = 0; i < narms; ++i){
      sumpos = 0;
      
      // loop over simulations
      for (int s = 0; s < nsims; ++s){
         nb = k[j] * ie[i * nsims + s] - ic[i * nsims + s];
        if (nb > 0.0){
          ++sumpos;
        }
      }
      prob.push_back((double)sumpos/nsims);
      k_vec.push_back(k[j]);
    }
  }
  return prob;
}


// Probability most effective therapy
//' @export
// [[Rcpp::export]]
std::vector<double> mceC(std::vector<double> k, std::vector<double> e,
                 std::vector<double> c,
                 int nsims, int narms) {
  
  int n_k = k.size();
  int N = n_k * narms;
  std::vector<double> prob(N, 0.0);
  
  // loop over willingess to pay k
  for (int j = 0; j < n_k; ++j){
    
    // loop over simulations
    for (int s = 0; s < nsims; ++s){
      std::vector<double> nb;
      nb.reserve(narms);
      
      // loop over treatment arms
      for (int i = 0; i < narms; ++i){
        nb.push_back(k[j] * e[s * narms + i] - c[s * narms + i]);
      }
      int nb_ind = vecmax_index(nb);
      prob[j * narms + nb_ind] = prob[j * narms + nb_ind] + 1;
    }
  }
  
  // Convert sum to proportion
  for (int i = 0; i < N; ++i){
    prob[i] = prob[i]/nsims;
  }
  
  return prob;
}

// Calculate Vstar - expected benefits if choosing optimal arm at each simulation
//' @export
// [[Rcpp::export]]
std::vector<double> VstarC(std::vector<double> k,
                          std::vector<double> e, std::vector<double> c,
                          int nsims, int narms) {

  int n_k = k.size();
  std::vector<double> Vstar;
  Vstar.reserve(n_k);

  // loop over willingess to pay k
  for (int j = 0; j < n_k; ++j){
    std::vector<double> Ustar;
    Ustar.reserve(nsims);
    
    // loop over simulations
    for (int s = 0; s < nsims; ++s){
      std::vector<double> ib;
      ib.reserve(narms);

      // loop over treatment arms
      for (int i = 0; i < narms; ++i){
        ib.push_back(k[j] * e[s * narms + i] - c[s * narms + i]);
      }
      Ustar.push_back(vecmax(ib));
    }
    double Ustar_mean = accumulate(Ustar.begin(), Ustar.end(), 0.0)/Ustar.size(); 
    Vstar.push_back(Ustar_mean);
  }

  return Vstar;
}

