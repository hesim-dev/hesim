// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/stats/survfuns.h>

// [[Rcpp::export]]
double C_test_rsurv(std::vector<double> time, std::vector<double> cumhaz,
                     bool time_inf = true){
  return hesim::stats::surv_sample(time, cumhaz, time_inf);
}