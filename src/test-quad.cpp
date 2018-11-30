#include <hesim/math/quad.h>
#include <hesim/check_R_infinity.h>


// [[Rcpp::export]]
double test_quad_dnorm(double lower, double upper){
  auto f = [](double x){
    return R::dnorm4(x, 0, 1, 0);
  };
  hesim::check_R_infinity(lower);
  hesim::check_R_infinity(upper);
  double abserr;
  int ier;
  return hesim::math::quad(f, lower, upper, abserr, ier);
}

// [[Rcpp::export]]
double test_quad_ier1(){
  auto f = [](double x){
    return 1/x;
  };
  double abserr;
  int ier;
  return hesim::math::quad(f, 1, INFINITY, abserr, ier);
}

// [[Rcpp::export]]
double test_quad_ier4(){
  auto f = [](double x){
    return sin(x);
  };
  double abserr;
  int ier;
  return hesim::math::quad(f, 1, INFINITY, abserr, ier);
}

// [[Rcpp::export]]
double test_quad_ier5(){
  auto f = [](double x){
    return 1/pow(x, 3);
  };
  double abserr;
  int ier;
  return hesim::math::quad(f, -2, 3, abserr, ier);
}


