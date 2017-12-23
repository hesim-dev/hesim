# ifndef DISTRIBUTIONS_H
# define DISTRIBUTIONS_H
#include <RcppArmadillo.h>
#include "utils.h"

vecmats convert_distribution_parameters(std::string dist, Rcpp::List R_parlist);
double qgompertzC (double p, double shape, double rate);
double rgompertzC (double shape, double rate);
double rsurv(double location, double anc1, std::string dist, double anc2 = 0.0);

class Distribution{
public:
  virtual ~Distribution() {};
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
  double pdf(double x) const;
  double cdf(double x) const;
  double quantile(double p) const;
  double hazard(double x) const;
  double cumhazard(double x) const;
  double random() const;
};

Distribution * select_distribution(std::string dist_name, 
                                   std::vector<double> parameters);
# endif


