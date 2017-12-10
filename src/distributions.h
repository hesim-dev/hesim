# ifndef DISTRIBUTIONS_H
# define DISTRIBUTIONS_H

double qgompertzC (double p, double shape, double rate);

double rgompertzC (double shape, double rate);

double rsurv(double location, double anc1, std::string dist, double anc2 = 0.0);

class Distribution{
public:
  virtual ~Distribution() {};
  virtual double pdf(double x) {return 0.0;}
  virtual double cdf(double x) {return 0.0;}
  virtual double quantile(double p){return 0.0;}
  virtual double hazard(double x) {return 0.0;}
  virtual double cumhazard(double x) {return 0.0;}
  virtual double random() {return 0.0;}
};

class Exponential : public Distribution {
private: 
  double rate_;
  
public:
  Exponential(double rate);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  double hazard(double x);
  double cumhazard(double x);
  double random();
};

class Weibull : public Distribution {
private:
  double shape_;
  double scale_;

public:
  Weibull(double shape, double scale);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  double hazard(double x);
  double cumhazard(double x);
  double random();
};

class Gamma : public Distribution {
private:
  double shape_;
  double rate_;
  
public:
  Gamma(double shape, double rate);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  double hazard(double x);
  double cumhazard(double x);
  double random();
};

class Lognormal : public Distribution {
private:
  double meanlog_;
  double sdlog_;
  
public:
  Lognormal(double meanlog, double sdlog);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  double hazard(double x);
  double cumhazard(double x);
  double random();
};

class Gompertz : public Distribution {
private:
  double shape_;
  double rate_;
  
public:
  Gompertz(double shape, double rate);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  double hazard(double x);
  double cumhazard(double x);
  double random();
};

class LogLogistic : public Distribution {
private:
  double shape_;
  double scale_;
  
public:
  LogLogistic(double shape, double scale);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  double hazard(double x);
  double cumhazard(double x);
  double random();
};

class GeneralizedGamma : public Distribution {
private:
  double mu_;
  double sigma_;
  double Q_;
  
public:
  GeneralizedGamma(double mu, double sigma, double Q);
  double pdf(double x);
  double cdf(double x);
  double quantile(double p);
  double hazard(double x);
  double cumhazard(double x);
  double random();
};

# endif


