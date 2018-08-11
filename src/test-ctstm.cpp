#include <hesim/ctstm/ctstm.h>


/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<bool> C_test_is_absorbing(arma::mat m){
  hesim::ctstm::trans_mat tmat(m);
  return tmat.absorbing_;
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<int> C_test_trans_mat_trans_id(arma::mat m, int from_state){
  hesim::ctstm::trans_mat tmat(m);
  return tmat.trans_id(from_state);
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<int> C_test_trans_mat_to(arma::mat m, int from_state){
  hesim::ctstm::trans_mat tmat(m);
  return tmat.to(from_state);
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
int C_test_trans_mat_n_trans(arma::mat m){
  hesim::ctstm::trans_mat tmat(m);
  return tmat.n_trans_;
}