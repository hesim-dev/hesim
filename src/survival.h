# ifndef SURVIVAL_H
# define SURVIVAL_H
#include "distributions.h"

class Survival {
private:
  std::string dist_name_;
  vecmats surv_coefs_;
  vecmats surv_X_;
  int nsims_; int nindivs_; int npars_;
  
public:
  Survival(std::string dist_name, Rcpp::List surv_coefs, Rcpp::List surv_X);
  int get_nsims();
  int get_nindivs();
  std::vector<double> predict_surv_pars(int sim, int id);
  template <typename Func>
  arma::mat main_loop(Func fun, std::vector<double> xvec);
  arma::mat quantiles(std::vector<double> q); 
  arma::mat hazard(std::vector<double> t); 
  arma::mat cumhazard(std::vector<double> t); 
  arma::mat surv(std::vector<double> t); 
};

class HazardFun {
public:
  double operator()(Distribution * dist, double t) const { 
    return dist->hazard(t); 
  }
};

class CumhazardFun {
public:
  double operator()(Distribution * dist, double t) const { 
    return dist->cumhazard(t); 
  }
};

class SurvFun {
public:
  double operator()(Distribution * dist, double t) const { 
    return 1 - dist->cdf(t); 
  }
};

class QuantileFun {
public:
  double operator()(Distribution * dist, double q) const { 
    return dist->quantile(q); 
  }
};

# endif


