#include <Rcpp.h>
using namespace Rcpp;


// Gompertz Distribution
// [[Rcpp::export]]
double qGompertz (double p, double shape, double rate) {
  double q = 0;
  if (shape == 0){
    q = R::qexp(p, rate, 1, 0);
  }
  else if (shape < 0){
    q = INFINITY;
  }
  else {
    q = 1/shape * log(1 - shape * log(1 - p)/rate);
  }
  return q;
}

// [[Rcpp::export]]
double rGompertz (double shape, double rate){
  double u = R::runif(0,1);
  return qGompertz(u, shape, rate);
}

// Random Survival Times
// [[Rcpp::export]]
double rSurv(double location, double par2, std::string dist) {
  double surv = 0.0;
  if (dist == "exp"){
    double rate = exp(location);
    surv = R::rexp(rate);
  }
  else if (dist == "weibull"){
    double shape = par2;
    double scale = exp(location);
    surv = R::rweibull(shape, scale);
  }
  else if (dist == "gompertz"){
    double shape = par2;
    double rate = exp(location);
    surv = rGompertz(shape, rate);
  }
  return surv;
}

/*** R
n <- 1000
r1 <- replicate(n, rGompertz(1, 1))
r2 <- rgompertz(1000, 1, 1)
summary(r1)
summary(r2)
  */


