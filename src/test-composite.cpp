# include <hesim/math/composite.h>
# include <vector>

// [[Rcpp::export]]
double C_test_trapz(std::vector<double> x, std::vector<double> y){
  return hesim::math::trapz(x.begin(), x.end(), y.begin());
}



