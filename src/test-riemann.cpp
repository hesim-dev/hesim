#include <hesim/math/riemann.h>
#include <cmath>

// [[Rcpp::export]]
double test_riemann_x2(std::vector<double> x){
  auto f = [](double x){
    return pow(x, 2);
  };
  return hesim::math::riemann(x.begin(), x.end(), f);
}

// [[Rcpp::export]]
std::vector<double> test_cum_riemann_x2(std::vector<double> x){
  auto f = [](double x){
    return pow(x, 2);
  };
  return hesim::math::cum_riemann(x.begin(), x.end(), f);
}








