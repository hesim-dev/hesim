# ifndef DISTRIBUTIONS_H
# define DISTRIBUTIONS_H
#include <RcppArmadillo.h>
#include "utils.h"
#include <RcppNumerical.h>
#include "zeroin.h"
#include <memory>

namespace hesim{

/********************
* Distribution class
********************/
class Distribution{
public:
  virtual ~Distribution() {};
  virtual void set_params(std::vector<double> params){};
  virtual double pdf(double x) const {return 0.0;}
  virtual double cdf(double x) const {return 0.0;}
  virtual double quantile(double p) const {return 0.0;}
  virtual double hazard(double x) const {return 0.0;}
  virtual double cumhazard(double x) const {return 0.0;}
  virtual double random() const {return 0.0;}
};

class Exponential : public Distribution {
private: 
  double rate_;
  
public:
  Exponential(double rate);
  void set_params(std::vector<double> params);
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class Weibull : public Distribution {
private:
  double shape_;
  double scale_;

public:
  Weibull(double shape, double scale);
  static Weibull create_from_Nma(double a0, double a1);
  void set_params(std::vector<double> params);
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class WeibullNma : public Distribution {
private:
  Weibull wei_;
  Weibull create_from_Nma(double a0, double a1);

public:
  WeibullNma(double a0, double a1);
  void set_params(std::vector<double> params);
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class Gamma : public Distribution {
private:
  double shape_;
  double rate_;
  
public:
  Gamma(double shape, double rate);
  void set_params(std::vector<double> params);
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class Lognormal : public Distribution {
private:
  double meanlog_;
  double sdlog_;
  
public:
  Lognormal(double meanlog, double sdlog);
  void set_params(std::vector<double> params);
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class Gompertz : public Distribution {
private:
  double shape_;
  double rate_;
  
public:
  Gompertz(double shape, double rate);
  void set_params(std::vector<double> params);
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class LogLogistic : public Distribution {
private:
  double shape_;
  double scale_;
  
public:
  LogLogistic(double shape, double scale);
  void set_params(std::vector<double> params);
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class GeneralizedGamma : public Distribution {
private:
  double mu_;
  double sigma_;
  double Q_;
  
public:
  GeneralizedGamma(double mu, double sigma, double Q);
  void set_params(std::vector<double> params);
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class SurvSplines : public Distribution {
private:
  std::vector<double> gamma_;
  std::vector<double> knots_;
  std::string scale_;
  std::string timescale_;
  int n_knots_;
  double knot_max_;
  double knot_min_;
  double timescale_fun(double x) const;
  double timescale_dx_fun(double x) const;
  double basis_cube(double x) const;
  double basis_cube_dx(double x) const;

public:
  SurvSplines(std::vector<double> gamma, std::vector<double> knots,
              std::string scale, std::string timescale);
  void set_params(std::vector<double> params);
  double linear_predict(double x) const;
  double linear_predict_dx(double x) const;
  double survival(double x) const;
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

class FracPoly : public Distribution {
private:
  std::vector<double> gamma_;
  std::vector<double> powers_;
  double basis_power(double x, double power) const;
  std::vector<double> basis(double x) const;
  
public:
  FracPoly(std::vector<double> gamma, std::vector<double> powers);
  void set_params(std::vector<double> params);
  double linear_predict(double x) const;
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

Distribution * select_distribution(std::string dist_name,
                                   std::vector<double> parameters);

std::unique_ptr<Distribution> init_Distribution_ptr(std::string dist_name,
                                   std::vector<double> parameters);


/*****************************
* Free functions and functors
*****************************/
class HazardFunc: public Numer::Func {
  private:
    const Distribution * dist_;
  public:
    HazardFunc(const Distribution * dist)
      : dist_(dist){}
    double operator()(const double& x) const {
      return dist_->hazard(x);
    }
};

inline double integrate_hazard(const Distribution * dist, double t){
    hesim::HazardFunc fun(dist);
    const double lower = 0, upper = t;
    double err_est; int err_code;
    return Numer::integrate(fun, lower, upper, err_est, err_code);
};

class ZeroinFunc {
  private:
    const Distribution * dist_;
    double p_;
  public:
    ZeroinFunc(const Distribution * dist, double p)
      : dist_(dist), p_(p){}
    double operator()(const double& x) const {
      return dist_->cdf(x) - p_;
    }
};

inline double quantile_numeric_work(const hesim::Distribution * dist, double p){
    hesim::ZeroinFunc func(dist, p);
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
    return zeroin(lower, upper, f_lower, f_upper, func,
                  &tol, &maxiter);
}

inline double quantile_numeric(const hesim::Distribution * dist, double p){
  if ( p < 0 || p > 1){
    return NAN;
  }
  else if (p == 0) return R_NegInf;
  else if (p == 1) return R_PosInf;
  else{
    return quantile_numeric_work(dist, p);
  }
}

class QuantileNumericFunc {
  private:
    const Distribution * dist_;
  public:
    QuantileNumericFunc(const Distribution * dist)
      : dist_(dist){}
    double operator()(const double& p) const {
      return quantile_numeric(dist_, p);
    }
};

class DiscountSurvFunc: public Numer::Func {
private:
  Distribution * dist_;
  double r_;
public:
  DiscountSurvFunc(Distribution * dist, double r) 
    : dist_(dist), r_(r){}
  
  double operator()(const double& x) const {
    return exp(-r_ * x) * (1 - dist_->cdf(x));
  }
};

inline double rmst(Distribution * dist, double t, double r){
  DiscountSurvFunc fun(dist, r);
  const double lower = 0, upper = t;
  double err_est; int err_code;
  return Numer::integrate(fun, lower, upper, err_est, err_code);
}

vecmats convert_distribution_parameters(std::string dist, Rcpp::List R_parlist);
double qgompertzC (double p, double shape, double rate);
double rgompertzC (double shape, double rate);
double rsurv(double location, double anc1, std::string dist, double anc2 = 0.0);

} //end namespace hesim

# endif


