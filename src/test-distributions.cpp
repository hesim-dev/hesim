#include <hesim/stats/distributions.h>
#include "test-distributions.h"


/**
 * \ingroup test
 * Function to test hesim::stats::exponential.
 */
// [[Rcpp::export]]
double C_test_exponential(double rate, std::string fun, double x = 0){
  hesim::stats::exponential exp(rate);
  return test_distribution(exp, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::weibull.
 */
// [[Rcpp::export]]
double C_test_weibull(double shape, double scale, std::string fun, double x = 0){
  hesim::stats::weibull wei(shape, scale);
  return test_distribution(wei, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::weibull_nma.
 */
// [[Rcpp::export]]
double C_test_weibull_nma(double a0, double a1, std::string fun, double x = 0){
  hesim::stats::weibull_nma wei(a0, a1);
  return test_distribution(wei, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::gamma.
 */
// [[Rcpp::export]]
double C_test_gamma(double shape, double rate, std::string fun, double x = 0){
  hesim::stats::gamma gamma(shape, rate);
  return test_distribution(gamma, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::lognormal.
 */
// [[Rcpp::export]]
double C_test_lognormal(double meanlog, double sdlog, std::string fun, double x = 0){
  hesim::stats::lognormal lnorm(meanlog, sdlog);
  return test_distribution(lnorm, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::gompertz.
 */
// [[Rcpp::export]]
double C_test_gompertz(double shape, double rate, std::string fun, double x = 0){
  hesim::stats::gompertz gompertz(shape, rate);
  return test_distribution(gompertz, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::loglogistic.
 */
// [[Rcpp::export]]
double C_test_loglogistic(double shape, double scale, std::string fun, double x = 0){
  hesim::stats::loglogistic llogis(shape, scale);
  return test_distribution(llogis, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::gengamma.
 */
// [[Rcpp::export]]
double C_test_gengamma(double mu, double sigma, double Q, std::string fun, double x = 0){
  hesim::stats::gengamma gengamma(mu, sigma, Q);
  return test_distribution(gengamma, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::survspline.
 */
// [[Rcpp::export]]
double C_test_survspline(std::vector<double> gamma, std::vector<double> knots, 
                         std::string scale, std::string timescale,
                         std::string fun, double x = 0){
  hesim::stats::survspline spline(gamma, knots, scale, timescale);
  return test_distribution(spline, fun, x);
}

/**
 * \ingroup test
 * Function to test hesim::stats::fracpoly.
 */
// [[Rcpp::export]]
double C_test_fracpoly(std::vector<double> gamma, std::vector<double> powers, 
                         std::string fun, double x = 0){
  hesim::stats::fracpoly fp(gamma, powers);
  return test_distribution(fp, fun, x);
}


/**
 * \ingroup test
 * Function to test hesim::stats::rtruncnorm.
 */
// [[Rcpp::export]]
double C_test_rtruncnorm(double mean, double sd, double lower, double upper){
  return hesim::stats::rtruncnorm(mean, sd, lower, upper);
}