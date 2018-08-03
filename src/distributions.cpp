#include <hesim/distributions.h>

/***************************************************************************//**
 * @ingroup stats
 * Vectorized random number generation for the generalized gamma distribution.
 * A vectorized version of hesim::stats::gengamma.random that is exported to @c R and
 * used in the @c R function @c fast_rgengamma.
 ******************************************************************************/ 
// [[Rcpp::export]]
std::vector<double> C_rgengamma(int n, std::vector<double> mu,
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
      hesim::stats::gengamma gengamma(mu[index], sigma[index], Q[index]);
      sample[i] = gengamma.random();
  }
  return sample;
}

/***************************************************************************//**
 * @ingroup stats
 * Vectorized random number generation for the piecewise exponential distribution.
 * A vectorized version of hesim::stats::rpwexp that is exported to @c R and
 * used in the @c R function @c rpwexp.
 ******************************************************************************/ 
// [[Rcpp::export]]
std::vector<double> C_rpwexp(int n, arma::mat rate, arma::rowvec time) {
  int b = rate.n_rows;
  std::vector<double> surv;
  surv.reserve(n);
  for (int i = 0; i < n; ++i){
      surv.push_back(hesim::stats::rpwexp(rate.row(i % b), time));
  }
  return surv;
}

/***************************************************************************//**
 * @ingroup stats
 * Vectorized random number generation for categorical distribution.
 * A vectorized version of hesim::stats::rcat that is exported to @c R and
 * used in the @c R function @c rcat.
 ******************************************************************************/ 
// [[Rcpp::export]]
std::vector<double> C_rcat(int n, arma::mat probs){
  int b = probs.n_rows;
  std::vector<double> samples(n);
  for (int i = 0; i < n; ++i){
    samples[i] = hesim::stats::rcat(probs.row(i % b));
  }
  return(samples);
}

/***************************************************************************//**
 * @ingroup stats
 * Vectorized random number generation for the Dirichlet distribution.
 * A vectorized version of hesim::stats::rdirichlet that is exported to @c R and
 * used in the @c R function @c rdirichlet_mat.
 ******************************************************************************/ 
// [[Rcpp::export]]
arma::cube C_rdirichlet_mat(int n, arma::mat alpha){
  int J = alpha.n_rows;
  int K = alpha.n_cols;
  arma::cube samp(J, K, n);
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < J; ++j){
      samp.slice(i).row(j) = hesim::stats::rdirichlet(alpha.row(j));
    }
  }
  return(samp);
}
/***************************************************************************//**
 * @ingroup stats
 * Function to test hesim::stats::rtruncnorm.
 ******************************************************************************/ 
// [[Rcpp::export]]
double C_test_rtruncnorm(double mean, double sd, double lower, double upper){
  return hesim::stats::rtruncnorm(mean, sd, lower, upper);
}

/***************************************************************************//**
 * @ingroup stats
 * Rcpp modules containing classes inherited from hesim::stats:::distribution.
 ******************************************************************************/ 
RCPP_MODULE(distributions){
  Rcpp::class_<hesim::stats::distribution>("distribution")
  .method("pdf", &hesim::stats::distribution::pdf)
  .method("cdf", &hesim::stats::distribution::cdf)
  .method("quantile", &hesim::stats::distribution::quantile)
  .method("hazard", &hesim::stats::distribution::hazard)
  .method("cumhazard", &hesim::stats::distribution::cumhazard)
  .method("random", &hesim::stats::distribution::random)
  ;

  Rcpp::class_<hesim::stats::exponential>("exponential")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double>()
    .method("pdf", &hesim::stats::exponential::pdf)
    .method("cdf", &hesim::stats::exponential::cdf)
    .method("quantile", &hesim::stats::exponential::quantile)
    .method("hazard", &hesim::stats::exponential::hazard)
    .method("cumhazard", &hesim::stats::exponential::cumhazard)
    .method("random", &hesim::stats::exponential::random)
  ;

  Rcpp::class_<hesim::stats::weibull>("weibull")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::stats::weibull::pdf)
    .method("cdf", &hesim::stats::weibull::cdf)
    .method("quantile", &hesim::stats::weibull::quantile)
    .method("hazard", &hesim::stats::weibull::hazard)
    .method("cumhazard", &hesim::stats::weibull::cumhazard)
    .method("random", &hesim::stats::weibull::random)
  ;
  
  Rcpp::class_<hesim::stats::weibull_nma>("weibull_nma")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::stats::weibull_nma::pdf)
    .method("cdf", &hesim::stats::weibull_nma::cdf)
    .method("quantile", &hesim::stats::weibull_nma::quantile)
    .method("hazard", &hesim::stats::weibull_nma::hazard)
    .method("cumhazard", &hesim::stats::weibull_nma::cumhazard)
    .method("random", &hesim::stats::weibull_nma::random)
  ;
  
  Rcpp::class_<hesim::stats::gamma>("gamma")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::stats::gamma::pdf)
    .method("cdf", &hesim::stats::gamma::cdf)
    .method("quantile", &hesim::stats::gamma::quantile)
    .method("hazard", &hesim::stats::gamma::hazard)
    .method("cumhazard", &hesim::stats::gamma::cumhazard)
    .method("random", &hesim::stats::gamma::random)
  ;
  
  Rcpp::class_<hesim::stats::lognormal>("lognormal")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::stats::lognormal::pdf)
    .method("cdf", &hesim::stats::lognormal::cdf)
    .method("quantile", &hesim::stats::lognormal::quantile)
    .method("hazard", &hesim::stats::lognormal::hazard)
    .method("cumhazard", &hesim::stats::lognormal::cumhazard)
    .method("random", &hesim::stats::lognormal::random)
  ;
  
  Rcpp::class_<hesim::stats::gompertz>("gompertz")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::stats::gompertz::pdf)
    .method("cdf", &hesim::stats::gompertz::cdf)
    .method("quantile", &hesim::stats::gompertz::quantile)
    .method("hazard", &hesim::stats::gompertz::hazard)
    .method("cumhazard", &hesim::stats::gompertz::cumhazard)
    .method("random", &hesim::stats::gompertz::random)
  ;
  
  Rcpp::class_<hesim::stats::loglogistic>("loglogistic")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::stats::loglogistic::pdf)
    .method("cdf", &hesim::stats::loglogistic::cdf)
    .method("quantile", &hesim::stats::loglogistic::quantile)
    .method("hazard", &hesim::stats::loglogistic::hazard)
    .method("cumhazard", &hesim::stats::loglogistic::cumhazard)
    .method("random", &hesim::stats::loglogistic::random)
  ;
  
  Rcpp::class_<hesim::stats::gengamma>("gengamma")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double, double, double>()
    .method("pdf", &hesim::stats::gengamma::pdf)
    .method("cdf", &hesim::stats::gengamma::cdf)
    .method("quantile", &hesim::stats::gengamma::quantile)
    .method("hazard", &hesim::stats::gengamma::hazard)
    .method("cumhazard", &hesim::stats::gengamma::cumhazard)
    .method("random", &hesim::stats::gengamma::random)
  ;
  
    Rcpp::class_<hesim::stats::survspline>("survspline")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<std::vector<double>, std::vector<double>, std::string, std::string>()
    .method("linear_predict", &hesim::stats::survspline::linear_predict)
    .method("linear_predict_dx", &hesim::stats::survspline::linear_predict_dx)
    .method("pdf", &hesim::stats::survspline::pdf)
    .method("cdf", &hesim::stats::survspline::cdf)
    .method("quantile", &hesim::stats::survspline::quantile)
    .method("hazard", &hesim::stats::survspline::hazard)
    .method("cumhazard", &hesim::stats::survspline::cumhazard)
    .method("random", &hesim::stats::survspline::random)
  ;
    
     Rcpp::class_<hesim::stats::fracpoly>("fracpoly")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<std::vector<double>, std::vector<double> >()
    .method("linear_predict", &hesim::stats::fracpoly::linear_predict)
    .method("pdf", &hesim::stats::fracpoly::pdf)
    .method("cdf", &hesim::stats::fracpoly::cdf)
    .method("quantile", &hesim::stats::fracpoly::quantile)
    .method("hazard", &hesim::stats::fracpoly::hazard)
    .method("cumhazard", &hesim::stats::fracpoly::cumhazard)
    .method("random", &hesim::stats::fracpoly::random)
  ;
}
