#include <hesim/quad.h>
#include <Rcpp.h>

class test_functor{
public:
  double operator()(double x){
    return 2* x;
  }
};

// [[Rcpp::export]]
double test_quad_functor(){
  test_functor f;
  double abserr;
  int ier;
  return hesim::math::quad(f, 0, 2, abserr, ier);
}

// [[Rcpp::export]]
double test_quad_lambda(){
  auto f = [](double x){
    return pow(x, 2);
  };
  double abserr;
  int ier;
  return hesim::math::quad(f, 0, 2, abserr, ier);
}
