// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/stats/survfuns.h>

// [[Rcpp::export]]
double C_test_rsurv(std::vector<double> time, std::vector<double> est,
                    std::string type = "surv", bool time_inf = true){
  return hesim::stats::surv_sample(time, est, type, time_inf);
}