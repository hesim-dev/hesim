# ifndef SURVIVAL_H
# define SURVIVAL_H
#include "distributions.h"
#include <RcppNumerical.h>

class Survival {
private:
  std::string dist_name_;
  vecmats surv_coefs_;
  vecmats surv_X_;
  int nsims_; int nindivs_; int npars_; int index_;
  double r_health_; double r_costs_;
  
public:
  Survival(std::string dist_name, Rcpp::List surv_coefs, Rcpp::List surv_X,
           double r_health, double r_costs);
  void set_index(int index);
  std::vector<double> predict_surv_pars(int sim, int id);
  void summary1(int sim, int id, std::string type, Distribution * dist, arma::mat &output, 
                std::vector<double> xvec);
  void main_loop(std::string type, std::vector<double> xvec, arma::mat &output);
  arma::mat summary(std::vector<double> xvec, std::string type);
};

class RmstFunc: public Numer::Func {
private:
  Distribution * dist_;
  double r_;
public:
  RmstFunc(Distribution * dist, double r) : dist_(dist), r_(r){}
  
  double operator()(const double& x) const {
    return exp(-r_ * x) * (1 - dist_->cdf(x));
  }
};

double rmst(Distribution * dist, double t, double r){
  RmstFunc fun(dist, r);
  const double lower = 0, upper = t;
  double err_est; int err_code;
  return Numer::integrate(fun, lower, upper, err_est, err_code);
}

// class QalysFunc: public Numer::Func {
// private:
//   double r_;
//   std::vector<double> w_;
//   Distribution * dist
// public:
//   QalysFun(Distribution * dist, double r, std::vector<double> w) : r_(r), w_(w) {}
//   
//   double operator()(const double& x) const {
//     return w_[x] * exp(-r_ * x) * (1 - dist->cdf(x));
// 
//   }
// };

class QalysFunc: public Numer::Func {
private:
  double r_;
  double w_;
  Exponential exponential_;
public:
  QalysFunc(double r, double w, Exponential exponential) : 
    r_(r), w_(w), exponential_(exponential){}

  double operator()(const double& x) const {
    return w_ * exp(-r_ * x) * (1 - exponential_.cdf(x));

  }
};

# endif


