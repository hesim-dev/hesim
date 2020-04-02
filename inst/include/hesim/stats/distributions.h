# ifndef DISTRIBUTIONS_H
# define DISTRIBUTIONS_H
#include <hesim/utils.h>
#include <hesim/Rbase/zeroin.h>
#include <hesim/math/quad.h>
#include <hesim/stats/survfuns.h>
#include <hesim/stats/rtrunc.h>
#include <memory>

namespace hesim{

/** @ingroup stats 
 * Probability distributions, random number 
 * generation, and statistical functions.
 */
namespace stats{


/***************************************************************************//** 
 * An abstract base class for probability distributions.
 * Each child probability distribution class contains a number of statistical 
 * functions such as a probability density function, 
 * cumulative distribution function, quantile function, hazard function, 
 * cumulative hazard function, and random number generator.
 ******************************************************************************/ 
class distribution{
public:
  virtual ~distribution() {};
  
  double max_x_ = INFINITY; ///< Maximum value for support of distribution. Can be infinity. 
                           ////< Only used when randomly sampling using numerical methods by 
                          ////< sampling from a cumulative hazard function.
  std::string cumhaz_method_; ///< Method used to compute the cumulative hazard 
                             ///< when is must be done numerically by integrating the hazard function.
  double step_; ////< Step size used to compute cumulative hazard at discrete points. 
               ////< Only used when randomly sampling using numerical methods by 
               ////< sampling from a cumulative hazard function.
  
  /** 
   * Set the parameters for the distribution.
   * @param param A vector of the parameters.
   * @return None.
   */ 
  virtual void set_params(std::vector<double> params) = 0;
  
  /** 
   * Probability density function.
   * @param x Quantile of the distribution.
   * @return. The density evaluated at @p x.
   */   
  virtual double pdf(double x) const = 0;

  /** 
   * Cumulative density function.
   * @param x Quantile of the distribution.
   * @return The distribution function evaluated at @p x.
   */     
  virtual double cdf(double x) const = 0;
  
  /** 
   * Quantile function.
   * @param p A probability to calculate a quantile for.
   * @return The quantile evaluated at @p p.
   */      
  virtual double quantile(double p) const = 0;
  
  /** 
   * Hazard function.
   * @param x Quantile of the distribution.
   * @return The hazard rate evaluated at @p x.
   */    
  virtual double hazard(double x) const {return 0.0;}
  
  /** 
   * Cumulative hazard function.
   * @param x Quantile of the distribution.
   * @return The cumulative hazard evaluated at @p x. 
   */   
  virtual double cumhazard(double x) const {return 0.0;}
  
  /** 
   * Random number generator.
   * @return A random sample from the probability distribution.
   */     
  virtual double random() const = 0;
  
  /**
   * Random number generator for a truncated distribution.
   * @return A random sample between a lower and upper bound from the probability distribution.
   */
  virtual double trandom(double lower, double upper) const {
    Rcpp::stop("The selected distribution cannot randomly sample from a truncated distribution.");
  }
  
}; // end base class distribution

/** Internal details for hesim::stats that should be ignored by external users.*/
namespace detail {

inline double quantile_numeric_work(const stats::distribution * dist, double p){
    auto func = [dist, p](const double& x){ // A lambda to pass to zeroin()
      return dist->cdf(x) - p;
    };
    double lower = -1;
    double upper = 1;
    while(func(lower) * func(upper) >= 0){
        double interval = upper - lower;
        lower = lower - 0.5 * interval;
        upper = upper + 0.5 * interval;
    }
    double f_lower = func(lower);
    double f_upper = func(upper);
    double tol = 0.0001;
    int maxiter = 1000;
    return Rbase::zeroin(lower, upper, f_lower, f_upper, func,
                  &tol, &maxiter);
};

inline double random_numeric(const stats::distribution * dist, 
                             std::string random_method) {
  if (random_method == "invcdf"){
    return dist->quantile(R::runif(0, 1));
  }
  else if (random_method == "discrete") {
    return surv_sample(dist, 0, INFINITY, dist->max_x_);
  }
  else{
    Rcpp::stop("'random_method' must be either 'invcdf' or 'discrete'.");
  }
}

inline double trandom_numeric(const stats::distribution * dist, 
                              double lower, double upper,
                              std::string random_method) {
  if (random_method == "invcdf"){
    return rtrunc(dist, lower, upper, "invcdf");
  }
  else if (random_method == "sample"){
    return rtrunc(dist, lower, upper, "sample"); 
  }
  else{
    Rcpp::stop("'random_method' must be either 'invcdf' or 'sample'.");
  }
}

} // end namespace detail


/***************************************************************************//** 
 * Compute quantile numerically.
 * Uses numerical methods to calculate the quantile of a probability distribution.
 * Should be used for distributions where quantiles cannot be computed 
 * analytically. Equivalent to the function @c qgeneric in the @c R package
 * @c flexsurv.
 * @param dist A pointer to the base class of a probability distribution.  
 * @param p A probability to calculate a quantile for.
 * @return The quantile evaluated at @p p. 
 ******************************************************************************/ 
inline double quantile_numeric(const distribution * dist, double p){
  if ( p < 0 || p > 1){
    return NAN;
  }
  else if (p == 0) return R_NegInf;
  else if (p == 1) return R_PosInf;
  else{
    return detail::quantile_numeric_work(dist, p);
  }
}


/***************************************************************************//**
 * The exponential distribution.
 * Uses the same parameterization as in the @c R package @c stats.
 ******************************************************************************/ 
class exponential : public distribution {
private: 
  double rate_; ///< rate parameter
  
public:
  /** 
   * The constructor.
   * Instantiates an exponential distribution with a given rate parameter.
   */ 
  exponential(double rate){
    rate_ = rate;
  }
  
  void set_params(std::vector<double> params){
    rate_ = exp(params[0]);
  }
  
  double pdf(double x) const{
    return rate_ * exp(-rate_ * x);
  }
  
  double cdf(double x) const{
    return 1 - exp(-rate_ * x); // R::pexp(x_, 1/rate_, 1, 0)
  }
  
  double quantile(double p) const{
    return R::qexp(p, 1/rate_, 1, 0);
  }
  
  double hazard(double x) const{
    return rate_;
  }
  
  double cumhazard(double x) const{
    return rate_ * x;
  }
  
  double random() const {
    return R::rexp(1/rate_);
  }
  
  double trandom(double lower, double upper) const   {
    return rtrunc(this, lower, upper, "invcdf");
  }
  
}; // end class exponential

/***************************************************************************//**
 * The Weibull distribution.
 * Uses the same parameterization as in the @c R package @c stats.
 ******************************************************************************/ 
class weibull : public distribution {
private:
  double shape_; ///< shape parameter
  double scale_; ///< scale parameter

public:
  /** 
   * The constructor.
   * Instantiates a Weibull distribution with a shape and scale parameters.
   */ 
  weibull(double shape, double scale) {
    shape_ = shape;
    scale_ = scale;
  }
  
  void set_params(std::vector<double> params) {
    shape_ = exp(params[0]);
    scale_ = exp(params[1]);
  }
  
  double pdf(double x) const {
    return R::dweibull(x, shape_, scale_, 0);
  }
  
  double cdf(double x) const {
    return R::pweibull(x, shape_, scale_, 1, 0);
  }
  
  double quantile(double p) const {
    return R::qweibull(p, shape_, scale_, 1, 0);
  }
  
  double hazard(double x) const {
    return shape_ * pow(x/scale_, shape_ - 1)/scale_;
  }
  
  double cumhazard(double x) const {
    return pow(x/scale_, shape_);
  }
  
  double random() const {
    return R::rweibull(shape_, scale_);
  }
  
  double trandom(double lower, double upper) const   {
    return rtrunc(this, lower, upper, "invcdf");
  }  
}; // end class weibull

/***************************************************************************//**
 * The Weibull distribution parameterized for network meta-analysis.
 ******************************************************************************/ 
class weibull_nma : public distribution {
private:
  weibull wei_; ///< A member of class @link weibull.

  /** 
   * Create an @link weibull object from a Weibull distribution parameterized
   * for a network meta-analysis.
   */ 
  weibull create_from_Nma(double a0, double a1){
    double shape = a1 + 1;
    double scalePH = exp(a0)/shape;
    double scale = pow(scalePH, -1/shape);
    return weibull(shape, scale);
  }

public:
  /** 
   * The constructor.
   * Instantiates a Weibull distribution with parameters @p a0 and @p a1.
   */  
  weibull_nma(double a0, double a1)
    : wei_(create_from_Nma(a0, a1)) {
  }
  
  void set_params(std::vector<double> params) {
    wei_ = create_from_Nma(params[0], params[1]);
  }
  
  double pdf(double x) const {
    return wei_.pdf(x);
  }
  
  double cdf(double x) const {
    return wei_.cdf(x);
  }
  
  double quantile(double p) const {
    return wei_.quantile(p);
  }
  
  double hazard(double x) const {
    return wei_.hazard(x);
  }
  
  double cumhazard(double x) const {
    return wei_.cumhazard(x);
  }
  
  double random() const {
    return wei_.random();
  }
  
  double trandom(double lower, double upper) const   {
    return rtrunc(this, lower, upper, "invcdf");
  }
  
}; // end class weibull_nma

/***************************************************************************//**
 * The gamma distribution.
 * Uses the same parameterization as in the @c R package @c stats..
 ******************************************************************************/ 
class gamma : public distribution {
private:
  double shape_; ///< Shape parameter.
  double rate_; ///< Rate parameter.
  
public:
  /** 
   * The constructor.
   * Instantiates a gamma distribution with a shape and rate parameters.
   */ 
  gamma(double shape, double rate) {
    shape_ = shape;
    rate_ = rate; 
  }
  
  void set_params(std::vector<double> params) {
    shape_ = exp(params[0]);
    rate_ = exp(params[1]);
  }
  
  double pdf(double x) const {
    return R::dgamma(x, shape_, 1/rate_, 0);
  }
  
  double cdf(double x) const {
    return R::pgamma(x, shape_, 1/rate_, 1, 0);
  }
  
  double quantile(double p) const {
    return R::qgamma(p, shape_, 1/rate_, 1, 0);
  }
  
  double hazard(double x) const {
    return gamma::pdf(x)/(1 - gamma::cdf(x));
  }
  
  double cumhazard(double x) const {
    return -R::pgamma(x, shape_, 1/rate_, 0, 1);
  }
  
  double random() const {
    return R::rgamma(shape_, 1/rate_);
  }
  
  double trandom(double lower, double upper) const   {
    return rtrunc(this, lower, upper, "invcdf");
  }
  
}; //end class gamma

/***************************************************************************//**
 * The lognormal distribution.
 * Uses the same parameterization as in the @c R package @c stats.
 ******************************************************************************/ 
class lognormal : public distribution {
private:
  double meanlog_; ///< Mean on log scale.
  double sdlog_; ///< Standard deviation on log scale.
  
public:
  /** 
   * The constructor.
   * Instantiates a lognormal distribution with mean and standard deviation parameters
   * defined on the log scale.
   */ 
  lognormal(double meanlog, double sdlog){
    meanlog_ = meanlog;
    sdlog_ = sdlog;
  }
  
  void set_params(std::vector<double> params) {
    meanlog_ = params[0];
    sdlog_ = exp(params[1]);
  }
  
  double pdf(double x) const {
    return R::dlnorm(x, meanlog_, sdlog_, 0);
  }
  
  double cdf(double x) const {
    return R::plnorm(x, meanlog_, sdlog_, 1, 0);
  }
  
  double quantile(double p) const {
    return R::qlnorm(p, meanlog_, sdlog_, 1, 0);
  }
  
  double hazard(double x) const {
    return lognormal::pdf(x)/(1 - lognormal::cdf(x));
  }
  
  double cumhazard(double x) const {
    return -R::plnorm(x, meanlog_, sdlog_, 0, 1);
  }
  
  double random() const {
    return R::rlnorm(meanlog_, sdlog_);
  }
  
  double trandom(double lower, double upper) const   {
    return rtrunc(this, lower, upper, "invcdf");
  }  
}; //end class lognormal

/***************************************************************************//**
 * The gompertz distribution.
 * Uses the same parameterization as in the @c R package @c flexsurv.
 ******************************************************************************/ 
class gompertz : public distribution {
private:
  double shape_; ///< Shape parameter.
  double rate_; ///< Rate parameter.
  
public:
  /** 
   * The constructor.
   * Instantiates a gompertz distribution with shape and rate parameters.
   */ 
  gompertz(double shape, double rate) {
    shape_ = shape;
    rate_ = rate;
  }
  
  void set_params(std::vector<double> params) {
    shape_ = params[0];
    rate_ = exp(params[1]); 
  }
  
  double pdf(double x) const {
    if (shape_ == 0){
      return R::dexp(x, 1/rate_, 0);
    }
    else{
      return rate_ * exp(shape_ * x) * exp(-rate_/shape_ * (exp(shape_ * x) -1));
    }
  }
  
  double cdf(double x) const {
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
  
  double quantile(double p) const {
    double asymp = 1 - exp(rate_/shape_);
    if (shape_ == 0){
        return R::qexp(p, 1/rate_, 1, 0);
    }
    else if (shape_ < 0 && p > asymp){
        return INFINITY;
    }
    else {
        return 1/shape_ * log(1 - shape_ * log(1 - p)/rate_);
    }
  }
  
  double hazard(double x) const {
    return rate_ * exp(shape_ * x);
  }
  
  double cumhazard(double x) const {
    if (shape_ == 0){
      return rate_ * x;
    }
    else{
      return rate_/shape_ * expm1(shape_ * x);
    }
  }
  
  double random() const {
    double u = R::runif(0,1);
    return quantile(u);
  }
  
  double trandom(double lower, double upper) const   {
    return rtrunc(this, lower, upper, "invcdf");
  }  
}; //end class gompertz

/***************************************************************************//**
 * The log-logistic distribution.
 * Uses the same parameterization as in the @c R package @c flexsurv.
 ******************************************************************************/ 
class loglogistic : public distribution {
private:
  double shape_; ///< Shape parameter.
  double scale_; ///< Scale parameter.
  
public:
  /** 
   * The constructor.
   * Instantiates a log-logistic distribution with shape and scale parameters.
   */ 
  loglogistic(double shape, double scale) {
    shape_ = shape;
    scale_ = scale;
  }
  
  void set_params(std::vector<double> params) {
    shape_ = exp(params[0]);
    scale_ = exp(params[1]);  
  }
  
  double pdf(double x) const {
    return (shape_/scale_) * pow((x/scale_), shape_ - 1)/
    pow((1 + pow((x/scale_), shape_)), 2); 
  }
  
  double cdf(double x) const {
    return 1 - 1/(1 + pow(x/scale_, shape_));
  }
  
  double quantile(double p) const {
    return exp(R::qlogis(p, log(scale_), 1/shape_, 1, 0)); 
  }
  
  double hazard(double x) const {
    return (shape_/scale_) * pow((x/scale_), shape_ - 1)/
    (1 + pow((x/scale_), shape_));
  }
  
  double cumhazard(double x) const {
    return -log(1 - loglogistic::cdf(x));
  }
  
  double random() const {
    double u = R::runif(0,1);
    return quantile(u);
  }
  
  double trandom(double lower, double upper) const   {
    return rtrunc(this, lower, upper, "invcdf");
  }  
}; //end class loglogistic


/***************************************************************************//**
 * The generalized gamma distribution.
 * Uses the same parameterization as in the @c R package @c flexsurv.
 ******************************************************************************/ 
class gengamma : public distribution {
private:
  double mu_; ///< Location parameter.
  double sigma_; ///< Scale parameter.
  double Q_; ///< Shape parameter.
  
public:
  /** 
   * The constructor.
   * Instantiates a generalized gamma distribution with location, shape, and scale
   * parameters.
   */ 
  gengamma(double mu, double sigma, double Q) {
    mu_ = mu;
    sigma_ = sigma;
    Q_ = Q;
  }
  
  void set_params(std::vector<double> params) {
    mu_ = params[0];
    sigma_ = exp(params[1]);
    Q_ = params[2];
  }
  
  double pdf(double x) const {
    if (Q_ != 0){
      double y = log(x);
      double w = (y - mu_)/sigma_;
      double Q2inv = 1/(Q_ * Q_);
      double logp = -log(sigma_ * x) + log(std::abs(Q_)) + Q2inv * log(Q2inv) +
      Q2inv * (Q_ * w - exp(Q_ * w)) - R::lgammafn(Q2inv);
      return exp(logp);
    } 
    else{
      return R::dlnorm(x, mu_, sigma_, 0);
    }
  }
  
  double cdf(double x) const {
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
  
  double quantile(double p) const {
    if (Q_ == 0){
      return R::qlnorm(p, mu_, sigma_, 1, 0);
    }
    else {
      if (Q_ < 0){
        p = 1 - p;
      }
      double gamma_quantile = R::qgamma(p, 1/(Q_ * Q_), 1, 1, 0);
      return exp(mu_ + sigma_ * (log(Q_ * Q_ * gamma_quantile)/Q_));
    }    
  }
  
  double hazard(double x) const {
    return gengamma::pdf(x)/(1 - gengamma::cdf(x));
  }
  
  double cumhazard(double x) const {
    return -log(1 - gengamma::cdf(x));
  }
  
  double random() const {
    if (Q_ == 0.0){
        return R::rlnorm(mu_, sigma_);
    }
    else{
        double w = log(pow(Q_, 2) * R::rgamma(1/pow(Q_, 2), 1))/Q_;
        return exp(mu_ + sigma_ * w);
    }
  }
  
  double trandom(double lower, double upper) const   {
    return rtrunc(this, lower, upper, "invcdf");
  }  
  
}; // end class gengamma

/***************************************************************************//**
 * A spline survival distribution.
 * Uses the same parameterization as in the @c R package @c flexsurv. Based on
 * Royston and Parmar (2002)
 ******************************************************************************/ 
class survspline : public distribution {
private:
  std::vector<double> gamma_; ///< The parameters of the spline model with length
                              ///< equal to the number of knots.
  std::vector<double> knots_; ///< Location of knots (including boundary knots).
  std::string scale_; ///< The transformation of the survival function determined by
                     ///< choice of g(). Options are "log_cumhazard", for the 
                     ///< log cummulative hazard; "log_hazard" for the log hazard 
                     ///< rate; "log_cumodds" for the log cummulative odds; and 
                     ///< "inv_normal" for the inverse normal distribution function.
  std::string timescale_; ///< If "log" (the default), then survival is modeled as a 
                          ///< spline function of log time; if "identity", then 
                          ///< it is modeled as a spline function of time.
  int n_knots_; ///< Number of knots.
  double knot_max_; ///< The largest knot.
  double knot_min_; ///< The smallest knot.
  std::string random_method_; ///< Method used to randomly draw from 
                             ///< survival function.                            
  
  // Function of time used for modeling survival.
  double timescale_fun(double x) const {
    if (timescale_ == "log"){
      return log(x);
    }
    else if (timescale_ == "identity"){
      return x;
    }
    else{
      Rcpp::stop("Selected timescale is not available.");
    }  
  }
  
  // Derivative of the function of time used for modeling survival.
  double timescale_dx_fun(double x) const {
    if (timescale_ == "log"){
        return 1/x;
    }
    else if (timescale_ == "identity"){
        return 1;
    }
    else{
        Rcpp::stop("Selected timescale is not available.");
    }
  }
  
  // A cubic basis function.
  double basis_cube(double x) const {
    if (x <= 0) {
      return 0;
    }
    else {
      return x * x * x;
    }    
  }
  
  // Derivative of a cubic basis function.
  double basis_cube_dx(double x) const {
    if (x <= 0) {
      return 0;
    }
    else {
      return 3 * x * x;
    }
  }

public:
  
  /** 
   * The constructor.
   * Instantiates a spline survival distribution.
   */ 
  survspline(std::vector<double> gamma, std::vector<double> knots,
              std::string scale, std::string timescale,
              std::string cumhaz_method = "quad", double step = -1,
              std::string random_method = "invcdf") {
    if (gamma.size() != knots.size()){
      Rcpp::stop("Length of gamma should equal number of knots.");
    }
    gamma_ = gamma;
    knots_ = knots;
    scale_ = scale;
    timescale_ = timescale;
    n_knots_ = knots.size();
    knot_max_ = *(knots.end() - 1);
    knot_min_ = *(knots.begin());
    cumhaz_method_ = cumhaz_method;
    step_ = step;
    random_method_ = random_method;
  }
  
  void set_params(std::vector<double> params) {
    gamma_ = params;
  }

   // Compute the spline function with a linear predictor (i.e., the spline function).
  double linear_predict(double x) const {
    double x_scaled = timescale_fun(x);
    std::vector<double> basis(n_knots_);
    basis[0] = 1; basis[1] = x_scaled;
    for (int j = 1; j < n_knots_ - 1; ++j){
      double lambda_j = (knot_max_ - knots_[j]) / (knot_max_ - knot_min_);
      basis[j + 1] = basis_cube(x_scaled - knots_[j]) - lambda_j *  basis_cube(x_scaled - knot_min_) -
      (1 - lambda_j) * basis_cube(x_scaled - knot_max_);
    }
    return std::inner_product(gamma_.begin(), gamma_.end(), basis.begin(), 0.0);
  }
  
  // Compute the derivative of the linear predictor (i.e., the spline function).
  double linear_predict_dx(double x) const {
    double x_scaled = timescale_fun(x);
    std::vector<double> basis_dx(n_knots_);
    basis_dx[0] = 0; basis_dx[1] = 1;
    for (int j = 1; j < n_knots_ - 1; ++j){
      double lambda_j = (knot_max_ - knots_[j]) / (knot_max_ - knot_min_);
      basis_dx[j + 1] = basis_cube_dx(x_scaled - knots_[j]) - lambda_j *  basis_cube_dx(x_scaled - knot_min_) -
      (1 - lambda_j) * basis_cube_dx(x_scaled - knot_max_);
    }
    return std::inner_product(gamma_.begin(), gamma_.end(), basis_dx.begin(), 0.0); 
  }
  
  /** 
   * Survival function.
   */   
  double survival(double x) const {
    if (x <= 0){
      return 1; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard" || scale_ == "log_cumhazard"){
      return exp(-cumhazard(x));
    }
    else if (scale_ == "log_cumodds"){
      return 1/(1 + exp(linear_predict(x)));
    }
    else if (scale_ == "inv_normal"){
      return R::pnorm(-linear_predict(x), 0, 1, 1, 0);
    }
    else {
      Rcpp::stop("Selected scale is not available.");
    }
  }
  
  double pdf(double x) const {
    if (x <= 0){
      return 0; // spline model is for time >= 0
    }
    double lp = linear_predict(x);
    double prob;
    if (scale_ == "log_hazard"){
      prob = survival(x) * hazard(x);
    }
    else if (scale_ == "log_cumhazard"){
      prob = timescale_dx_fun(x) * linear_predict_dx(x) * exp(lp - exp(lp));
    }
    else if (scale_ == "log_cumodds"){
      prob = timescale_dx_fun(x) * linear_predict_dx(x) * exp(lp - 2 * log(1 + exp(lp)));
    }
    else if (scale_ == "inv_normal"){
        prob = timescale_dx_fun(x) * linear_predict_dx(x) * R::dnorm(lp, 0, 1, 0);
    }
    else{
      Rcpp::stop("Selected scale is not available.");
    }
    if (prob <= 0){
      prob = 0;
    }
    return prob;    
  }
  
  double cdf(double x) const {
    return 1 - survival(x);
  }
  
  double quantile(double p) const {
    return quantile_numeric(this, p);
  }
  
  double hazard(double x) const {
    if (x <= 0){
      return 0; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard"){
      return exp(linear_predict(x));
    }
    else if (scale_ == "log_cumhazard"){
      return timescale_dx_fun(x) * linear_predict_dx(x) * exp(linear_predict(x));
    }
    else if (scale_ == "log_cumodds"){
      return timescale_dx_fun(x) * linear_predict_dx(x) * R::plogis(linear_predict(x), 0, 1, 1, 0);
    }
    else if (scale_ == "inv_normal"){
      double lp = linear_predict(x);
      return timescale_dx_fun(x) * linear_predict_dx(x) * R::dnorm(-lp, 0, 1, 0)/R::pnorm(-lp, 0, 1, 1, 0);
    }
    else{
      Rcpp::stop("Selected scale is not available.");
    }
  }
  double cumhazard(double x) const {
    if (x <= 0){
      return 0; // spline model is for time >= 0
    }
    if (scale_ == "log_hazard"){
      return integrate_hazard(this, x, cumhaz_method_);
    }
    else if (scale_ == "log_cumhazard"){
      return exp(linear_predict(x));
    }
    else if (scale_ == "log_cumodds"){
      return log1p(exp(linear_predict(x)));
    }
    else if (scale_ == "inv_normal"){
      return -R::pnorm(-linear_predict(x), 0, 1, 1, 1);
    }
    else{
      Rcpp::stop("Selected scale is not available.");
    }
  }
  
  double random() const {
    return detail::random_numeric(this, random_method_);
  }
  
  double trandom(double lower, double upper) const   {
    return detail::trandom_numeric(this, lower, upper, random_method_);
  }    
};

/***************************************************************************//**
 * A fractional polynomial survival distribution.
 * Details to be provided...
 ******************************************************************************/ 
class fracpoly : public distribution {
private:
  std::vector<double> gamma_; ///< The scale and shape parameters.
  std::vector<double> powers_; ///< The powers of the fractional polynomial.
  std::string random_method_; ///< Method used to randomly draw from 
                             ///< survival function.                               
  
  // 
  double basis_power(double x, double power) const {
    if (power == 0){
      return log(x);
    }
    else{
      return pow(x, power);
    }    
  }
  
  std::vector<double> basis(double x) const {
    int n_powers = powers_.size();
    std::vector<double> basis(n_powers + 1);
    basis[0] = 1;
    basis[1] = basis_power(x, powers_[0]);
    double xp_old = basis[1];
    double xp_new;
    if (n_powers > 1){
      for (int i = 1; i < n_powers; ++i){
        if (powers_[i] == powers_[i - 1]){
          xp_new = log(x) * xp_old;
        }
        else {
          xp_new = basis_power(x, powers_[i]);
        }
        basis[i + 1] = xp_new;
        xp_old = xp_new;
      }
    }
    return basis;
  }
  
public:
   
  /** 
   * The constructor.
   * Instantiates a fractional polynomial survival distribution.
   */ 
  fracpoly(std::vector<double> gamma, std::vector<double> powers,
           std::string cumhaz_method = "quad", double step = -1,
           std::string random_method = "invcdf") {
    gamma_ = gamma;
    powers_ = powers;
    cumhaz_method_ = cumhaz_method;
    step_ = step;
    random_method_ = random_method;
  }
  
  void set_params(std::vector<double> params) {
    gamma_ = params;
  }
  
  double linear_predict(double x) const {
    std::vector<double> b = basis(x);
    return std::inner_product(gamma_.begin(), gamma_.end(), b.begin(), 0.0);
  }
  
  double pdf(double x) const {
    return hazard(x) * (1 - cdf(x));  
  }
  
  double cdf(double x) const {
    return 1 - exp(-cumhazard(x));
  }
  
  double quantile(double p) const {
    return quantile_numeric(this, p); 
  }
  
  double hazard(double x) const {
    if (x <= 0){
      return 0; //  model is for time >= 0
    }
    else{
      return exp(linear_predict(x));
    }
  }
  
  double cumhazard(double x) const {
    return integrate_hazard(this, x, cumhaz_method_);
  }
  
  double random() const {
    return detail::random_numeric(this, random_method_);
  }
  
  double trandom(double lower, double upper) const   {
    return detail::trandom_numeric(this, lower, upper, random_method_);
  }    
};

/***************************************************************************//**
 * Random number generation for the categorical distribution.
 * @param probs A vector of length @c K specifying the probability of each of the 
 * @c K categories. Internally normalized to sum to 1.
 * @return A random sample from the categorical distribution.
 ******************************************************************************/ 
inline int rcat(arma::rowvec probs) {
  int k = probs.n_elem;
  double probs_sum = accu(probs);
  probs = probs/probs_sum;
  Rcpp::IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  int max = which_max(ans);
  return(max);
}

/***************************************************************************//**
 * Random number generation for the truncated normal distribution.
 * @param mean Mean of the distribution.
 * @param sd Standard deviation of the distribution.
 * @param lower Lower bound.
 * @param Upper bound
 * @return A random sample from the truncated normal distribution.
 ******************************************************************************/ 
inline double rtruncnorm(double mean, double sd, double lower, double upper){
  double  sample;
  sample = R::rnorm(mean, sd);
  while(sample < lower || sample > upper){
      sample = R::rnorm(mean, sd);
  }
  return sample;
}

/***************************************************************************//**
 * Random number generation for the piecewise exponential distribution.
 * @param rate Vector of rates with each element denoting a unique time
 * period.
 * @param time Vector of times of the same length as @p rate giving the times
 * at which the rate changes.
 * @return A random sample from the piecewise exponential distribution.
 ******************************************************************************/ 
inline double rpwexp (arma::rowvec rate, arma::rowvec time) {
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

/***************************************************************************//**
 * Random number generation for the Dirichlet distribution.
 * @param alpha Vector of concentration parameters.
 * @return A random sample from the Dirichlet distribution.
 ******************************************************************************/ 
inline arma::rowvec rdirichlet(arma::rowvec alpha){
  int alpha_len = alpha.size();
  arma::rowvec x(alpha_len);
  for (int i = 0; i < alpha_len; ++i){
      x(i) = R::rgamma(alpha(i), 1);
  }
  return x/arma::sum(x);
}

} //end namespace stats

} //end namespace hesim

# endif


