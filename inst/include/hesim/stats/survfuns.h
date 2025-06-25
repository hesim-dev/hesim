# ifndef HESIM_STATS_SURVFUNS_H
# define HESIM_STATS_SURVFUNS_H

#include <hesim/math/quad.h>
#include <hesim/math/riemann.h>
#include <hesim/utils.h>

namespace hesim{

namespace stats{

namespace detail {

/***************************************************************************//** 
 * Integrate a hazard function using quadrature.
 * Integrate a hazard function from 0 to t using Gaussian quadrature. Also see
 * quad().
 * @param dist A pointer to the base class of a probability distribution.  
 * @param t Time to integrate hazard until.
 * @return The integral of the hazard function.
 ******************************************************************************/ 
template <class Dist>
inline double integrate_hazard_quad(Dist dist, double t){
  auto fun = [dist](double x){
    return dist->hazard(x);
  };
  const double lower = 0, upper = t;
  double abserr; int ier;
  return math::quad(fun, lower, upper, abserr, ier);
};

/***************************************************************************//** 
 * Integrate a hazard function using a riemann sum.
 * Integrate a hazard function from 0 to t using a riemann sum. Also see
 * riemann().
 * @param dist A pointer to the base class of a probability distribution.  
 * @param t Time to integrate hazard until.
 * @return The integral of the hazard function.
 ******************************************************************************/ 
template <class Dist>
inline double integrate_hazard_riemann(Dist dist, double t){
  auto fun = [dist](double x){
    return dist->hazard(x);
  };
  if (t <= 0){
    return 0;
  }
  else{
    std::vector<double> times = hesim::seq(0, t, dist->step_);
    return math::riemann(times.begin(), times.end(), fun); 
  }
};

}

/***************************************************************************//**
 * Random number generation for the bernoulli distribution.
 * @param n Number of samples.
 * @param p Probability of success.
 * @return A random sample from the bernoulli distribution.
 ******************************************************************************/ 
inline int rbernoulli(double p){
  return R::runif(0, 1) > (1 - p);
}

/***************************************************************************//** 
 * Integrate a hazard function.
 * Integrate a hazard function from 0 to t using a variety of integration 
 * methods.
 * @param dist A pointer to the base class of a probability distribution.  
 * @param t Time to integrate hazard until.
 * @param method Integration method to use. Options are "quad" for quadrature
 * or "riemann" for an approximation via a riemann sum.
 * @return The integral of the hazard function.
 ******************************************************************************/ 
template <class Dist>
inline double integrate_hazard(Dist dist, double t, std::string method){
  if (method == "quad"){
    return detail::integrate_hazard_quad(dist, t);
  }
  else if (method == "riemann"){
    return detail::integrate_hazard_riemann(dist, t);
  }
  else {
    Rcpp::stop("The integration method must be 'quad' or 'riemann'.");
  }
};

/***************************************************************************//** 
 * Compute cumulative hazard numerically.
 * Compute a cumulative hazard function at discrete time points numerically by
 *  integrating the hazard function. 
 * @param hazfun A hazard function.  
 * @param times Times to compute the cumulative hazard function at. 
 * @param method Integration method to use. Options are "quad" for quadrature
 * or "riemann" for an approximation via a riemann sum.
 * @return The cumulative hazard.
 ******************************************************************************/ 
template <class Func>
inline std::vector<double> cumhazard_numeric(Func hazfun, std::vector<double> times,
                                             std::string method){
  if (method == "quad"){
    std::vector<double> cumhazard(times.size());
    const double lower = 0;
    double abserr; int ier;
    for (int i = 0; i < (int) times.size(); ++i){
      const double upper = times[i];
      cumhazard[i] = math::quad(hazfun, lower, upper, abserr, ier);
    }
    return cumhazard;  
  }
  else if (method == "riemann"){
    return math::cum_riemann(times.begin(), times.end(), hazfun);
  }
  else {
    Rcpp::stop("The integration method must be 'quad' or 'riemann'.");
  }
};

/***************************************************************************//** 
 * Compute restricted mean survival time.
 * Compute restricted mean survival time over a given time period for a chosen
 * probability distribution. Optionally discount survival at a rate @p r > 0. 
 * @param dist A pointer to the base class of a probability distribution.  
 * @param t Time to calculate mean survival time until.
 * @return Restricted mean survival time.
 ******************************************************************************/ 
template <class Dist>
inline double rmst(Dist dist, double t, double r = 0){
  auto fun = [dist, r](double x){
    return exp(-r * x) * (1 - dist->cdf(x));
  };
  const double lower = 0, upper = t;
  double err_est; int err_code;
  return math::quad(fun, lower, upper, err_est, err_code);
}

/***************************************************************************//** 
 * Randomly number generation from an arbitrary survival distribution.
 * Randomly draw a single observation from a survival distribution given 
 * discrete cumulative hazard curves or survival curves.
 * @param time Times at which estimates were computed. 
 * @param cumhaz Estimates of the cumulative hazard.
 * @param time_inf Determines whether the survival time of infinity be simulated. 
 * If true, then the probability of surviving to time INFINITY is assumed to equal 
 * the probability of surviving beyond the final time period in @p time; otherwise,
 * individuals are assumed to only survive to the final time period in @p time.
 * @return A random sample from the survival distribution.
 ******************************************************************************/ 
inline double surv_sample(std::vector<double> &time, std::vector<double> cumhaz,
                          bool time_inf = true){
  
  double died = 0;
  int i = 1;
  unsigned int n_times = time.size();
  while(died == 0 && i < (int) n_times){
    double prob = 1 - exp(cumhaz[i - 1] - cumhaz[i]);
    died = rbernoulli(prob);
    if (died == 1){
      return time[i];
    } 
    else{
      ++i; 
    }
  }
  return INFINITY;
}

/***************************************************************************//** 
 * Randomly number generation from an arbitrary survival distribution.
 * Randomly draw a single observation from a survival distribution given 
 * a hazard function. The hazard function is used to generate cumulative
 * hazard curves.
 * @param dist A pointer to the base class of a probability distribution.
 * @param lower, upper The cumulative hazard function is computed from @p lower 
 * to @p upper. @p lower must be non-negative.
 * @param max_survtime The maximum value of time that survival probabilities are
 * computed until, which must only be specified if @p upper equals INFINITY. In
 * this case, the cumulative hazard is computed until @p max_survtime, and the
 * probability of INFINITY is assumed to equal the probability of surviving beyond
 * @p max_survtime. Must be positive and cannot be infinite.
 * @return A random sample from the survival distribution.
 ******************************************************************************/ 
template <class Dist>
inline double surv_sample(Dist dist, double lower = 0, double upper = INFINITY,
                    double max_survtime = -1) {
  // Exceptions
  if (lower < 0){
    Rcpp::stop("'lower' cannot be negative.");
  }
  if (std::isinf(upper) && max_survtime < 0){
    Rcpp::stop("'max_survtime' cannot be negative.");
  }
  if (std::isinf(max_survtime)){
    Rcpp::stop("'max_survtime' cannot be infinite.");
  }
  
  // Times to compute hazards at
  std::vector<double> time;
  if (std::isinf(upper)){
    time = seq(lower, max_survtime, dist->step_);
    // double step = (1.0/12.0) * (max_survtime/100);
    // time = seq(lower, max_survtime, step);
  } 
  else{
    time = seq(lower, upper, dist->step_);
    // double step = (1.0/12.0) * (upper/100);
    // time = seq(lower, upper, step);
  }
  
  // Compute hazards
  auto hazfun = [dist](double x){
    return dist->hazard(x);
  };
  std::vector<double> cumhazard = cumhazard_numeric(hazfun, time, dist->cumhaz_method_);
  
  // Sample
  bool time_inf = false;
  if (std::isinf(upper)){
    time_inf = true;
  }
  return surv_sample(time, cumhazard, time_inf);
}

} // End namespace stats

} // End namespace hesim

# endif
