#include <hesim/ctstm/ctstm.h>

/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<bool> C_ctstm_is_absorbing(arma::mat m){
  return hesim::ctstm::is_absorbing(m);
}