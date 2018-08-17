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
