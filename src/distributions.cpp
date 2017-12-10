// [[Rcpp::interfaces(r, cpp)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "distributions.h"
using namespace Rcpp;

/**************************
* Exponential distribution
**************************/
Exponential::Exponential(double rate){
  rate_ = rate;
}

double Exponential::pdf(double x){
  return rate_ * exp(-rate_ * x);
}

double Exponential::cdf(double x){
  return 1 - exp(-rate_ * x); // R::pexp(x_, 1/rate_, 1, 0)
}

double Exponential::quantile(double p){
  return R::qexp(p, 1/rate_, 1, 0);
}

double Exponential::hazard(double x){
  return rate_;
}

double Exponential::cumhazard(double x){
  return rate_ * x;
}

double Exponential::random(){
  return R::rexp(1/rate_);
}

/*********************
* Weibull distribution
*********************/
Weibull::Weibull(double shape, double scale){
  shape_ = shape;
  scale_ = scale;
}

double Weibull::pdf(double x){
  return R::dweibull(x, shape_, scale_, 0);
}

double Weibull::cdf(double x){
  return R::pweibull(x, shape_, scale_, 1, 0);
}

double Weibull::quantile(double p){
  return R::qweibull(p, shape_, scale_, 1, 0);
}

double Weibull::hazard(double x){
  return shape_ * pow(x/scale_, shape_ - 1)/scale_;
}

double Weibull::cumhazard(double x){
  return pow(x/scale_, shape_);
}

double Weibull::random(){
  return R::rweibull(shape_, scale_);
}

/*******************
* Gamma distribution
*******************/
Gamma::Gamma(double shape, double rate){
  shape_ = shape;
  rate_ = rate;
}

double Gamma::pdf(double x){
  return R::dgamma(x, shape_, 1/rate_, 0);
}

double Gamma::cdf(double x){
  return R::pgamma(x, shape_, 1/rate_, 1, 0);
}

double Gamma::quantile(double p){
  return R::qgamma(p, shape_, 1/rate_, 1, 0);
}

double Gamma::hazard(double x){
  return Gamma::pdf(x)/(1 - Gamma::cdf(x));
}

double Gamma::cumhazard(double x){
  return -R::pgamma(x, shape_, 1/rate_, 0, 1);
}

double Gamma::random(){
  return R::rgamma(shape_, 1/rate_);
}

/***********************
* Lognormal distribution
***********************/
Lognormal::Lognormal(double meanlog, double sdlog){
  meanlog_ = meanlog;
  sdlog_ = sdlog;
}

double Lognormal::pdf(double x){
  return R::dlnorm(x, meanlog_, sdlog_, 0);
}

double Lognormal::cdf(double x){
  return R::plnorm(x, meanlog_, sdlog_, 1, 0);
}

double Lognormal::quantile(double p){
  return R::qlnorm(p, meanlog_, sdlog_, 1, 0);
}

double Lognormal::hazard(double x){
  return Lognormal::pdf(x)/(1 - Lognormal::cdf(x));
}

double Lognormal::cumhazard(double x){
  return -R::plnorm(x, meanlog_, sdlog_, 0, 1);
}

double Lognormal::random(){
  return R::rlnorm(meanlog_, sdlog_);
}

/**********************
* Gompertz distribution
**********************/
// [[Rcpp::export]]
double qgompertz(double p, double shape, double rate) {
  double asymp = 1 - exp(rate/shape);
  if (shape == 0){
    return R::qexp(p, 1/rate, 1, 0);
  }
  else if (shape < 0 && p > asymp){
    return INFINITY;
  }
  else {
    return 1/shape * log(1 - shape * log(1 - p)/rate);
  }
}

// [[Rcpp::export]]
double rgompertz(double shape, double rate){
  double u = R::runif(0,1);
  return qgompertz(u, shape, rate);
}

Gompertz::Gompertz(double shape, double rate){
  shape_ = shape;
  rate_ = rate;
}

double Gompertz::pdf(double x){
  if (shape_ == 0){
    return R::dexp(x, 1/rate_, 0);
  }
  else{
    return rate_ * exp(shape_ * x) * exp(-rate_/shape_ * (exp(shape_ * x) -1));
  }
}

double Gompertz::cdf(double x){
  if (shape_ == 0){
    return R::pexp(x, 1/rate_, 1, 0);
  }
  else if (std::isinf(x)){
    return 1;
  }
  else{
    return 1 - exp(-rate_/shape_ * (exp(shape_ * x) - 1));
  }
}

double Gompertz::quantile(double p){
  return qgompertz(p, shape_, rate_);
}

double Gompertz::hazard(double x){
  return rate_ * exp(shape_ * x);
}

double Gompertz::cumhazard(double x){
  if (shape_ == 0){
    return rate_ * x;
  }
  else{
    return rate_/shape_ * expm1(shape_ * x);
  }
}

double Gompertz::random(){
  return rgompertz(shape_, rate_);
}

/**************************
* Log-logistic distribution
**************************/
// [[Rcpp::export]]
double qllogis(double p, double shape, double scale, int lt = 1, int lg = 0){
  return exp(R::qlogis(p, log(scale), 1/shape, lt, lg));
}

// [[Rcpp::export]]
double rllogis(double shape, double scale){
  double u = R::runif(0,1);
  return qllogis(u, shape, scale);
}

LogLogistic::LogLogistic(double shape, double scale){
  shape_ = shape;
  scale_ = scale;
}

double LogLogistic::pdf(double x){
  return (shape_/scale_) * pow((x/scale_), shape_ - 1)/
    pow((1 + pow((x/scale_), shape_)), 2);
}

double LogLogistic::cdf(double x){
  return 1 - 1/(1 + pow(x/scale_, shape_));
}

double LogLogistic::quantile(double p){
  return exp(R::qlogis(p, log(scale_), 1/shape_, 1, 0));
}

double LogLogistic::hazard(double x){
  return (shape_/scale_) * pow((x/scale_), shape_ - 1)/
    (1 + pow((x/scale_), shape_));
}

double LogLogistic::cumhazard(double x){
  return -log(1 - LogLogistic::cdf(x));
}

double LogLogistic::random(){
  return rllogis(shape_, scale_);
}

/********************************
* Generalized gamma distribution
********************************/
// [[Rcpp::export]]
double rgengamma(double mu, double sigma, double Q){
  if (Q == 0.0){
    return R::rlnorm(mu, sigma);
  }
  else{
    double w = log(pow(Q, 2) * R::rgamma(1/pow(Q, 2), 1))/Q;
    return exp(mu + sigma * w);
  }
}

GeneralizedGamma::GeneralizedGamma(double mu, double sigma, double Q){
  mu_ = mu;
  sigma_ = sigma;
  Q_ = Q;
}

double GeneralizedGamma::pdf(double x){
  if (Q_ != 0){
    double y = log(x);
    double w = (y - mu_)/sigma_;
    double Q2inv = 1/(Q_ * Q_);
    double logp = -log(sigma_ * x) + log(std::abs(Q_)) + Q2inv * log(Q2inv) +
      Q2inv * (Q_ * w - exp(Q_ * w)) - R::lgammafn(Q2inv);
    return exp(logp);
  } // 
  else{
    return R::dlnorm(x, mu_, sigma_, 0);
  }
}

double GeneralizedGamma::cdf(double x){
  double y = log(x);
  double w = (y - mu_)/sigma_;
  double Q2inv = 1/(Q_ * Q_);
  double expnu = exp(Q_ * w) * Q2inv;
  if (Q_ > 0){
    return R::pgamma(expnu, Q2inv, 1, 1, 0);
  }
  else if (Q_ == 0){
    return R::plnorm(x, mu_, sigma_, 1, 0);
  }
  else{
    return 1 - R::pgamma(expnu, Q2inv, 1, 1, 0);
  }
}

double GeneralizedGamma::quantile(double p){
  if (Q_ == 0){
    return R::qlnorm(p, mu_, 1/(sigma_ * sigma_), 1, 0);
  }
  else {
    double gamma_quantile = R::qgamma(p, 1/(Q_ * Q_), 1, 1, 0);
    return exp(mu_ + sigma_ * (log(Q_ * Q_ * gamma_quantile)/Q_));
  }
}

double GeneralizedGamma::hazard(double x){
  return GeneralizedGamma::pdf(x)/(1 - GeneralizedGamma::cdf(x));
}

double GeneralizedGamma::cumhazard(double x){
  return -log(1 - GeneralizedGamma::cdf(x));
}

double GeneralizedGamma::random(){
  return rgengamma(mu_, sigma_, Q_);
}

// [[Rcpp::export(name="C_rgengamma_vec")]]
std::vector<double> rgengamma_vec(int n, std::vector<double> mu, 
                                  std::vector<double> sigma, 
                                  std::vector<double> Q){
  std::vector<double> sample(n);
  int mu_size = mu.size();
  int sigma_size = sigma.size();
  int Q_size = Q.size();
  if (mu_size != sigma_size || mu_size != Q_size){
    Rcpp::stop("Length of mu, sigma, and Q must be the same");
  }
  for (int i = 0; i < n; ++i){
    int index = i % mu_size;
    GeneralizedGamma gengamma(mu[index], sigma[index], Q[index]);
    sample[i] = gengamma.random(); 
  }
  return(sample);
}

/*******************************
* Truncated normal distribution
********************************/
// [[Rcpp::export]]
double rtruncnorm(double mean, double sd, double lower, double upper){
  double  sample;
  sample = R::rnorm(mean, sd);
  while(sample < lower || sample > upper){
    sample = R::rnorm(mean, sd);
  } 
  return sample;
}

/************************************
* Piecewise exponential distribution
*************************************/
// NOTE: rate in R::rexp is 1/rate in rexp!!!!!!!!
// [[Rcpp::export]]
double rpwexp (arma::rowvec rate, arma::rowvec time) {
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
// [[Rcpp::export(name="C_rpwexp_vec")]]
std::vector<double> rpwexp_vec (int n, arma::mat rate, arma::rowvec time) {
  int b = rate.n_rows;
  std::vector<double> surv;
  surv.reserve(n);
  for (int i = 0; i < n; ++i){
    surv.push_back(rpwexp(rate.row(i % b), time));
  }
  return surv;
}

/*************************
* Categorical distribution
**************************/
// [[Rcpp::export]]
int rcat(arma::rowvec probs) {
  int k = probs.n_elem;
  double probs_sum = accu(probs);
  probs = probs/probs_sum;
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  int max = which_max(ans);
  return(max);
}

// [[Rcpp::export(name="C_rcat_vec")]]
arma::vec rcat_vec(int n, arma::mat probs){
  int b = probs.n_rows;
  arma::vec samp(n);
  for (int i = 0; i < n; ++i){
    samp(i) = rcat(probs.row(i % b));
  }
  return(samp);
}

/***********************
* Dirichlet distribution
************************/
// [[Rcpp::export]]
arma::rowvec rdirichlet(arma::rowvec alpha){
  int alpha_len = alpha.size();
  arma::rowvec x(alpha_len);
  for (int i = 0; i < alpha_len; ++i){
    x(i) = R::rgamma(alpha(i), 1);
  }
  return x/arma::sum(x);
}

// [[Rcpp::export(name="C_rdirichlet_mat")]]
arma::cube rdirichlet_mat(int n, arma::mat alpha){
  int J = alpha.n_rows;
  int K = alpha.n_cols;
  arma::cube samp(J, K, n);
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < J; ++j){
      samp.slice(i).row(j) = rdirichlet(alpha.row(j));
    }
  }
  return(samp);
}
