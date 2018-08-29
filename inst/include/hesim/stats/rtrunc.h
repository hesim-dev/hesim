# ifndef HESIM_STATS_RTRUNC_H
# define HESIM_STATS_RTRUNC_H

# include <Rcpp/Rmath.h>

namespace hesim {
  
namespace stats {

/***************************************************************************//**
 * Randomly sample from a truncated probability distribution with repeated
 * sampling.
 * Generate random numbers by repeatedly sampling from the non-truncated 
 * distribution until the sample lies between the desired interval. Note that
 * this can be very slow if samples are nearly always drawn outside of the
 * interval.
 * @param f A functor or lambda expression used to randomly draw samples from
 * the non-truncated distribution.
 * @param lower, upper Lower and upper bounds of the random variable.
 * @return A random sample from the truncated normal distribution.
 ******************************************************************************/ 
template <typename Func>
inline double rtrunc_repeat(Func f, double lower, double upper){
  double sample = f();
  while(sample < lower || sample > upper){
      sample = f();
  }
  return sample;
}

/***************************************************************************//**
 * Compute quantile from a truncated probability distribution.
 * @param f_cdf A functor or lambda to compute the cumulative densitiy function 
 * of the non-truncated distribution.
 * @param f_quantile A functor or lambda to compute the quantile of the non-truncated
 * distribution.
 * @param p A probability to calculate a quantile for.
 * @param lower, upper Lower and upper bounds of the random variable.
 * @return The quantile evaluated at @p p.
 ******************************************************************************/ 
template <typename Func1, typename Func2>
inline double qtrunc(Func1 f_cdf, Func2 f_quantile, double p, double lower, double upper){
  if (f_cdf(lower) == f_cdf(upper)) {
    Rcpp::stop( "Truncation interval is not inside the domain of the quantile function");
  }    
  double v = f_cdf(lower) + (f_cdf(upper) - f_cdf(lower)) * p;
  return f_quantile(v);
}

/***************************************************************************//**
 * Randomly sample from a truncated probability distribution with inverse CDF.
 * Generate random numbers by using the inverse CDF method. Note that
 * this can fail if samples are required too far into the tail.
 * @param f A functor or lambda expression used to compute quantiles for the
 * the non-truncated distribution.
 * @param lower, upper Lower and upper bounds of the random variable.
 * @return A random sample from the truncated normal distribution.
 ******************************************************************************/ 
template <typename Func1, typename Func2>
inline double rtrunc_invcdf(Func1 f_cdf, Func2 f_quantile, double lower, double upper){
  double u = R::runif(0, 1);
  return qtrunc(f_cdf, f_quantile, u, lower, upper);
}

/***************************************************************************//**
 * Randomly sample from a truncated probability distribution.
 * Generate random numbers using using one of the methods available in @p method.
 * @param dist A probability distribution class. Must have member functions 
 * @c cdf for computing the cumulative density function, @c quantile for 
 * computing quantiles, and @c random for random number generation.
 * @param lower, upper Lower and upper bounds of the random variable.
 * @param method "invcdf" for the inverse CDF method as in hesim::stats::rtrunc_invcdf
 * and "repeat" for repeated sampling from the non-truncated distribution as in
 * hesim::stats::rtrunc_repeat. 
 * @return A random sample from the truncated normal distribution.
 ******************************************************************************/ 
template <class Dist>
inline double rtrunc(Dist dist, double lower, double upper, 
                     std::string method = "invcdf"){
  if (method == "invcdf"){
    auto f_cdf = [dist](double x){ return dist->cdf(x); };
    auto f_quantile = [dist](double p){ return dist->quantile(p); };    
    return rtrunc_invcdf(f_cdf, f_quantile, lower, upper);
  }
  else {
    auto f_random = [dist](){ return dist->random(); };
    return rtrunc_repeat(f_random, lower, upper);
  }
}


} // end namespace stats

} // end namespace hesim


# endif