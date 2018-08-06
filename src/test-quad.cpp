#include <hesim/math/quad.h>

class test_functor{
public:
  double operator()(double x){
    return R::dnorm4(x, 0, 1, 0);
  }
};

void check_infinity(double &x){
  if (!R_FINITE(x)){
    x = INFINITY;
  }
}

// [[Rcpp::export]]
double test_quad_functor(double lower, double upper){
  test_functor f;
  check_infinity(lower);
  check_infinity(upper);
  double abserr;
  int ier;
  return hesim::math::quad(f, lower, upper, abserr, ier);
}

// [[Rcpp::export]]
double test_quad_lambda(double lower, double upper){
  auto f = [](double x){
    return R::dnorm4(x, 0, 1, 0);
  };
  check_infinity(lower);
  check_infinity(upper);
  double abserr;
  int ier;
  return hesim::math::quad(f, lower, upper, abserr, ier);
}

