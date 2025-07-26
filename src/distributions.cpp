#include <hesim/stats/distributions.h>

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
std::vector<double> C_rpwexp(int n, arma::mat rate, std::vector<double> time) {
  int b = rate.n_rows;
  std::vector<double> out;
  out.reserve(n);
  
  for (int i = 0; i < n; ++i){
    arma::rowvec rate_i = rate.row(i % b);
    out.push_back(hesim::stats::rpwexp(rate_i, time));
  }
  return out;
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
 * Rcpp modules containing classes inherited from hesim::stats:::distribution.
 ******************************************************************************/ 
RCPP_MODULE(distributions){
  Rcpp::class_<hesim::stats::distribution>("distribution")
  .field("max_x_", &hesim::stats::distribution::max_x_)
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
    .method("trandom", &hesim::stats::exponential::trandom)
  ;
  
  Rcpp::class_<hesim::stats::piecewise_exponential>("piecewise_exponential")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<std::vector<double>, std::vector<double>>()
    .method("pdf", &hesim::stats::piecewise_exponential::pdf)
    .method("cdf", &hesim::stats::piecewise_exponential::cdf)
    .method("quantile", &hesim::stats::piecewise_exponential::quantile)
    .method("hazard", &hesim::stats::piecewise_exponential::hazard)
    .method("cumhazard", &hesim::stats::piecewise_exponential::cumhazard)
    .method("random", &hesim::stats::piecewise_exponential::random)
    .method("trandom", &hesim::stats::piecewise_exponential::trandom)
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
    .method("trandom", &hesim::stats::weibull::trandom)
  ;
  
  Rcpp::class_<hesim::stats::weibull_ph>("weibull_ph")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<double, double>()
    .method("pdf", &hesim::stats::weibull_ph::pdf)
    .method("cdf", &hesim::stats::weibull_ph::cdf)
    .method("quantile", &hesim::stats::weibull_ph::quantile)
    .method("hazard", &hesim::stats::weibull_ph::hazard)
    .method("cumhazard", &hesim::stats::weibull_ph::cumhazard)
    .method("random", &hesim::stats::weibull_ph::random)
    .method("trandom", &hesim::stats::weibull_ph::trandom)
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
    .method("trandom", &hesim::stats::weibull_nma::trandom)
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
    .method("trandom", &hesim::stats::gamma::trandom)
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
    .method("trandom", &hesim::stats::lognormal::trandom)
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
    .method("trandom", &hesim::stats::gompertz::trandom)
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
    .method("trandom", &hesim::stats::loglogistic::trandom)
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
    .method("trandom", &hesim::stats::gengamma::trandom)
  ;
  
    Rcpp::class_<hesim::stats::survspline>("survspline")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<std::vector<double>, std::vector<double>, std::string, 
                 std::string, std::string, double, std::string>()
    .method("linear_predict", &hesim::stats::survspline::linear_predict)
    .method("linear_predict_dx", &hesim::stats::survspline::linear_predict_dx)
    .method("pdf", &hesim::stats::survspline::pdf)
    .method("cdf", &hesim::stats::survspline::cdf)
    .method("quantile", &hesim::stats::survspline::quantile)
    .method("hazard", &hesim::stats::survspline::hazard)
    .method("cumhazard", &hesim::stats::survspline::cumhazard)
    .method("random", &hesim::stats::survspline::random)
    .method("trandom", &hesim::stats::survspline::trandom)
  ;
    
     Rcpp::class_<hesim::stats::fracpoly>("fracpoly")
    .derives<hesim::stats::distribution>("distribution")
    .constructor<std::vector<double>, std::vector<double>, std::string, 
                  double, std::string >()
    .method("linear_predict", &hesim::stats::fracpoly::linear_predict)
    .method("pdf", &hesim::stats::fracpoly::pdf)
    .method("cdf", &hesim::stats::fracpoly::cdf)
    .method("quantile", &hesim::stats::fracpoly::quantile)
    .method("hazard", &hesim::stats::fracpoly::hazard)
    .method("cumhazard", &hesim::stats::fracpoly::cumhazard)
    .method("random", &hesim::stats::fracpoly::random)
    .method("trandom", &hesim::stats::fracpoly::trandom)
  ;

   Rcpp::class_<hesim::stats::point_mass>("point_mass")
     .derives<hesim::stats::distribution>("distribution")
     .constructor<double>()
     .method("pdf", &hesim::stats::point_mass::pdf)
     .method("cdf", &hesim::stats::point_mass::cdf)
     .method("random", &hesim::stats::point_mass::random)
     .method("trandom", &hesim::stats::point_mass::trandom)
   ;  
}
