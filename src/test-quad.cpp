#include <hesim/math/quad.h>
#include <hesim/check_R_infinity.h>

class test_functor{
public:
  double operator()(double x){
    return R::dnorm4(x, 0, 1, 0);
  }
};

// [[Rcpp::export]]
double test_quad_functor(double lower, double upper){
  test_functor f;
  hesim::check_R_infinity(lower);
  hesim::check_R_infinity(upper);
  double abserr;
  int ier;
  return hesim::math::quad(f, lower, upper, abserr, ier);
}

// [[Rcpp::export]]
double test_quad_lambda(double lower, double upper){
  auto f = [](double x){
    return R::dnorm4(x, 0, 1, 0);
  };
  hesim::check_R_infinity(lower);
  hesim::check_R_infinity(upper);
  double abserr;
  int ier;
  return hesim::math::quad(f, lower, upper, abserr, ier);
}

