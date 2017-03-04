#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// Truncated normal distribution
// [[Rcpp::export]]
double rtruncnormC(double mean, double sd, double lower, double upper){
  double  sample;
  sample = R::rnorm(mean, sd);
  while(sample < lower || sample > upper){
    sample = R::rnorm(mean, sd);
  } 
  return sample;
}

// Gompertz distribution
// [[Rcpp::export]]
double qgompertzC(double p, double shape, double rate) {
  double q = 0;
  double asymp = 1 - exp(rate/shape);
  if (shape == 0){
    q = R::qexp(p, rate, 1, 0);
  }
  else if (shape < 0 && p > asymp){
    q = INFINITY;
  }
  else {
    q = 1/shape * log(1 - shape * log(1 - p)/rate);
  }
  return q;
}

// [[Rcpp::export]]
double rgompertzC(double shape, double rate){
  double u = R::runif(0,1);
  return qgompertzC(u, shape, rate);
}

// Log-logistic distribution
// [[Rcpp::export]]
double qllogisC(double p, double shape, double scale, int lt = 1, int lg = 0){
  return exp(R::qlogis(p, log(scale), 1/shape, lt, lg));
}

// [[Rcpp::export]]
double rllogisC(double shape, double scale){
  double u = R::runif(0,1);
  return qllogisC(u, shape, scale);
}

// Generalized gamma distribution
// [[Rcpp::export]]
double rgengammaC(double mu, double sigma, double Q){
  double samp = 0.0;
  if (Q == 0.0){
    samp = R::rlnorm(mu, sigma);
  }
  else{
    double w = log(pow(Q, 2) * R::rgamma(1/pow(Q, 2), 1))/Q;
    samp = exp(mu + sigma * w);
  }
  return samp;
}

// Piecewise exponential  
// NOTE: rate in R::rexp is 1/rate in rexp!!!!!!!!
// [[Rcpp::export]]
double rpwexp1C (arma::rowvec rate, arma::rowvec time) {
  int T = rate.n_elem;
  double surv = 0.0;
  for (int t = 0; t < T; ++t){
    double rexp_t = R::rexp(1/rate(t));
    surv = time(t) + rexp_t;
    if (t < (T - 1)){
      if (surv < time(t + 1)){
        break;
      }
    }
  }
  return surv;
}

// Vectorized piecewise exponential 
// [[Rcpp::export]]
std::vector<double> rpwexpC (arma::mat rate, arma::rowvec time) {
  int N = rate.n_rows;
  std::vector<double> surv;
  surv.reserve(N);
  for (int i = 0; i < N; ++i){
    surv.push_back(rpwexp1C(rate.row(i), time));
  }
  return surv;
}

// Random survival times
// [[Rcpp::export]]
double rsurv(double location, double anc1, std::string dist, double anc2 = 0.0) {
  double surv = 0.0;
  if (dist == "exponential"){
    double rate = exp(location);
    surv = R::rexp(1/rate);
  }
  else if (dist == "weibull"){
    double shape = exp(anc1);
    double scale = exp(location);
    surv = R::rweibull(shape, scale);
  }
  else if (dist == "gompertz"){
    double shape = anc1;
    double rate = exp(location);
    surv = rgompertzC(shape, rate);
  }
  else if (dist == "lnorm"){
    double meanlog = location;
    double sdlog = exp(anc1);
    surv = R::rlnorm(meanlog, sdlog);
  }
  else if (dist == "gamma"){
    double rate = exp(location);
    double shape = exp(anc1);
    surv = R::rgamma(shape, 1/rate);
  }
  else if (dist == "llogis"){
    double scale = exp(location);
    double shape = exp(anc1);
    surv = rllogisC(shape, scale);
  }
  else if (dist == "gengamma"){
    double mu = location;
    double sigma = exp(anc1);
    double Q = anc2;
    surv = rgengammaC(mu, sigma, Q);
  }
  return surv;
}

// Random generation for categorical distribution
// [[Rcpp::export]]
int rcat1C(arma::rowvec probs) {
  int k = probs.n_elem;
  double probs_sum = accu(probs);
  probs = probs/probs_sum;
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  int max = which_max(ans);
  return(max);
}

// [[Rcpp::export]]
arma::vec rcatC(arma::mat probs){
  int n = probs.n_rows;
  arma::vec samp(n);
  for (int i = 0; i < n; ++i){
    samp(i) = rcat1C(probs.row(i));
  }
  return(samp);
}

