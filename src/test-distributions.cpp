
#include <hesim/stats/distributions.h>

/**
 * \ingroup test
 * Function to test hesim::stats::rtruncnorm.
 */
// [[Rcpp::export]]
double C_test_rtruncnorm(double mean, double sd, double lower, double upper){
  return hesim::stats::rtruncnorm(mean, sd, lower, upper);
}