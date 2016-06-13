#include <Rcpp.h>
using namespace Rcpp;

// Factorial for integer n
// [[Rcpp::export]]
unsigned int factorialC(unsigned int n)
{
  unsigned int ret = 1;
  for (unsigned int i = 1; i <= n; ++i)
    ret *= i;
  return ret;
}

// Continuous present value for constant a, i.e. integral of a * exp(-rt)
// [[Rcpp::export]]
double pv1(double t1, double t2, double a, double r){
  double pv;
  if (r == 0){
    pv = a * (t2 - t1);
  }
  else {
    pv = a/r * (exp(-r * t1) - exp(-r * t2));
  }
  return pv;
}

// Integral for present value of exponent, i.e. (t-h)^p * exp(-rt)
// [[Rcpp::export]]
double exponent_int(double r, unsigned int p, double t, double h) {
  double sum = 0;
  double integral = 0;
  if (r > 0){
    for (unsigned int i = 0; i <=p; ++i){
      sum += 1/pow(r, i + 1) * factorialC(p)/factorialC(p - i) * pow(t - h, p - i);
    }
    integral = -exp(-r * t) * sum;
  }
  else {
    integral = pow(t - h, p + 1)/(p+1);
  }
  return integral;
}

// Integral of continously discounted polynomial of degree p, i.e.
// [a*(t-h)^p + .. a*(t-h) + a)] * exp(-rt)
// [[Rcpp::export]]
double poly_int(double r, unsigned int p, double t, double h, std::vector<double> a){
  double sum = 0;
  for (unsigned int i = 0; i <=p; ++i){
    sum += a[i] * exponent_int(r, i, t, h);
  }
  return sum;
}

// Definite integral of continously discounted polynomial of degree p, i.e.present
// value of integral [a*(t-h)^p + .. a*(t-h) + a)] * exp(-rt) from t1 to t2
// [[Rcpp::export]]
double pv_poly(double r, unsigned int p, double t1,
               double t2, double h, std::vector<double> a){
  double upper = poly_int(r, p, t2, h, a);
  double lower = poly_int(r, p, t1, h, a);
  return upper - lower;
}


/*** R
r <- 0.03
h <- 10
avec <- c(1, 1.2, 1.4, 10)
  a <- t(matrix(avec))
  pvPoly(r, 3, 10, 20, h, a)
  t <- seq(10, 20, .01)
  library("pracma")
  y <- a[1] * exp(-r * t) + a[2] * (t - h) * exp(-r * t) + a[3] * (t-h)^2 * exp(-r * t) +
    a[4] * (t - h)^3 * exp(-r * t)
  trapz(t, y)

  n <- 1e6
a <- matrix(avec, nrow = n, ncol = length(avec), byrow=TRUE)
  t1 <- rep(10, n)
  t2 <- rep(20, n)
  h <- rep(h, n)
  ptm <- proc.time()
  tmp <- pvPoly(r, 3, t1, t2, h, a)
  head(tmp)
  proc.time() - ptm
  */
