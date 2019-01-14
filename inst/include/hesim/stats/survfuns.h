# ifndef HESIM_STATS_SURVFUNS_H
# define HESIM_STATS_SURVFUNS_H

#include <hesim/Rbase/sample.h>
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
  std::vector<double> times = hesim::seq(0, t, dist->step_);
  return math::riemann(times.begin(), times.end(), fun);
};

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
 * Compute cumulative hazards by taking the cumulative sum of the
 * hazard within each discrete time interval.
 * @param hazfun A hazard function.
 * @param t Time to compute the cumulative hazard until.
 * @return The cumulative hazard.
 ******************************************************************************/ 
template <class Func>
inline std::vector<double> cumhazard_discrete(Func hazfun, std::vector<double> time){
  std::vector<double> cumhazard(time.size());
  cumhazard[0] = 0;
  for (int i = 1; i < time.size(); ++i){
    double step = time[i] - time[i - 1];
    cumhazard[i] = step * hazfun(time[i]) + cumhazard[i - 1];
  }
  return cumhazard;
}

/***************************************************************************//** 
 * Randomly number generation from an arbitrary survival distribution.
 * Randomly draw a single observation from a survival distribution given 
 * cumulative hazard curves or survival curves.
 * @param time Times at which estimates were computed. 
 * @param est Estimates of the cumulative hazard or survival curves.
 * @param type Is the estimate a cumulative hazard curve (@c "cumhazard") or
 * a survival curve (@c "surv").
 * @param time_inf Determines whether the survival time of infinity be simulated. 
 * If true, then the probability of surviving to time INFINITY is assumed to equal 
 * the probability of surviving beyond the final time period in @p time; otherwise,
 * individuals are assumed to only survive to the final time period in @p time.
 * @return A random sample from the survival distribution.
 ******************************************************************************/ 
inline double rsurv(std::vector<double> &time, std::vector<double> est,
                   std::string type = "cumhazard", bool time_inf = true){

  auto diff = [](std::vector<double> x){
    std::vector<double> x_diff(x.size());
    x_diff[0] = 0;
    for (int i = 1; i < x.size(); ++i){
      x_diff[i] = x[i] - x[i - 1];
    }
    return x_diff;
  };
  
  std::vector<double> cdf(est.size());
  if (type == "cumhazard"){
    for (int i = 0; i < est.size(); ++i) cdf[i] = 1 - exp(-est[i]);
  } 
  else if (type == "surv"){
    for (int i = 0; i < est.size(); ++i) cdf[i] = 1 - est[i];
  }
  else{
    Rcpp::stop("'type' must either be 'cumhazard' or 'surv'.");
  }
  
  arma::vec prob = diff(cdf);
  if (time_inf){
    int prob_size = prob.size();
    prob.resize(prob_size + 1);
    prob(prob_size) = 1 - cdf[cdf.size() - 1]; // Survival probability at last time point
    time.push_back(INFINITY); 
  }
  int size = 1;
  return Rbase::sample(time, size, false, prob)[0];
}

/***************************************************************************//** 
 * Randomly number generation from an arbitrary survival distribution.
 * Randomly draw a single observation from a survival distribution given 
 * a hazard function. The hazard function is used to generate cumulative
 * hazard curves.
 * @param hazfun A functor or lambda expression to compute the hazard as a 
 * function of time.
 * @param lower, upper The cumulative hazard function is computed from @p lower 
 * to @p upper. @p lower must be non-negative.
 * @param max_survtime The maximum value of time that survival probabilities are
 * computed until, which must only be specified if @p upper equals INFINITY. In
 * this case, the cumulative hazard is computed until @p max_survtime, and the
 * probability of INFINITY is assumed to equal the probability of surviving beyond
 * @p max_survtime. Must be positive and cannot be infinite.
 * @return A random sample from the survival distribution.
 ******************************************************************************/ 
template <class Func>
inline double rsurv(Func hazfun, double lower = 0, double upper = INFINITY,
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
    double step = (1.0/12.0) * (max_survtime/100);
    time = seq(lower, max_survtime, step);
  } 
  else{
    double step = (1.0/12.0) * (upper/100);
    time = seq(lower, upper, step);
  }
  
  // Compute hazards
  std::vector<double> cumhazard = cumhazard_discrete(hazfun, time);
  
  // Sample
  bool time_inf = false;
  if (std::isinf(upper)){
    time_inf = true;
  }
  return rsurv(time, cumhazard, "cumhazard", time_inf);
}

} // End namespace stats

} // End namespace hesim

# endif