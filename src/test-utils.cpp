// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <hesim/utils.h>

/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<int> C_test_add_constant_int(std::vector<int> v, double value){
  hesim::add_constant(v, value);
  return v;
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<double> C_test_add_constant_double(std::vector<double> v, double value){
  hesim::add_constant(v, value);
  return v;
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
double C_test_pv(double z, double r, double t1, double t2){
  return(hesim::pv(z, r, t1, t2));
}

// [[Rcpp::export]]
std::vector<double> C_test_seq(double from, double to, double by){
  return hesim::seq(from, to , by);
}

// [[Rcpp::export]]
double C_test_max_lt(std::vector<double> v, double value){
  auto it = hesim::max_lt(v.begin(), v.end(), value);
  return *it;
}

// [[Rcpp::export]]
double C_test_find_interval(double x, std::vector<double> vec){
  return hesim::find_interval(x, vec);
}

// [[Rcpp::export]]
bool C_test_is_member_of(double x, std::vector<double> lookup){
  return hesim::is_member_of(x, lookup);
}
