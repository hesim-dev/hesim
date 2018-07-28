#include <hesim/integrate.h>

class IntegrateTestFunc{
  public:
    double operator()(double x) const {
      return x * x;
    }
};

// [[Rcpp::export]]
double C_test_trapzfun(std::vector<double> x){
  IntegrateTestFunc fun;
  return hesim::math::trapzfun(fun, x.begin(), x.end());
}

// [[Rcpp::export]]
double C_test_trapz(std::vector<double> x, std::vector<double> y){
  return hesim::math::trapz(x.begin(), x.end(), y.begin());
}

// [[Rcpp::export]]
std::vector<double> C_test_cumtrapzfun(std::vector<double> t){
  IntegrateTestFunc fun;
  return hesim::math::cumtrapzfun(fun, t.begin(), t.end());
}

// [[Rcpp::export]]
double C_test_simpsfun(std::vector<double> t){
  IntegrateTestFunc fun;
  return hesim::math::simpsfun(fun, t.begin(), t.end());
}

// [[Rcpp::export]]
std::vector<double> C_test_cumsimpsfun(std::vector<double> t){
  IntegrateTestFunc fun;
  return hesim::math::cumsimpsfun(fun, t.begin(), t.end());
}
