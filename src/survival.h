# ifndef SURVIVAL_H
# define SURVIVAL_H
#include "distributions.h"
#include <RcppNumerical.h>

class RmstFunc: public Numer::Func {
private:
  Distribution * dist_;
  double r_;
  double time_length_;
public:
  RmstFunc(Distribution * dist, double r, double time_length) 
    : dist_(dist), r_(r), time_length_(time_length){}
  
  double operator()(const double& x) const {
    return exp(-r_ * x * time_length_) * (1 - dist_->cdf(x));
  }
};

double rmst(Distribution * dist, double t, double r, double time_length){
  RmstFunc fun(dist, r, time_length);
  const double lower = 0, upper = t;
  double err_est; int err_code;
  return Numer::integrate(fun, lower, upper, err_est, err_code);
}

class StateValues {
public:
  virtual ~StateValues() {};
  StateValues() 
    : sim_(0), id_(0), component_(0){}
  void set_indices(int sim, int id, int component = 0){
    sim_ = sim;
    id_ = id;
    component_ = component;
  }
  virtual double operator()() const = 0;

protected:
  int sim_;
  int id_;
  int component_;
};

class UtilityDefault : public StateValues {
public:
  UtilityDefault(SEXP weights) 
    : StateValues(){
    Rcpp::Environment weights_env = Rcpp::as<Rcpp::Environment>(weights);
    weights_ = Rcpp::as<std::vector<double> >(weights_env["mean"]);
  }
  double operator()() const {
    return weights_[sim_];
  }
  
protected:
  std::vector<double> weights_;
};

class CostsDefault : public StateValues {
public:
  CostsDefault(SEXP costs) 
    : StateValues(){
    Rcpp::Environment costs_env = Rcpp::as<Rcpp::Environment>(costs);
    Rcpp::List costs_list = Rcpp::as<Rcpp::List>(costs_env["mean"]);
    costs_ = list_to_vecmats(costs_list);
    
  }
  double operator()() const {
    return costs_[component_](sim_, 0);
  }
  
protected:
  vecmats costs_;
};

class WeightedRmstFunc: public Numer::Func {
private:
  Distribution * dist_;
  double r_;
  double time_length_;
  StateValues * state_values_;
public:
  WeightedRmstFunc(Distribution * dist, double r, double time_length, StateValues * state_values) 
    : dist_(dist), r_(r), time_length_(time_length), state_values_(state_values){}
  
  double operator()(const double& x) const {
    return (*state_values_)() * exp(-r_ * x * time_length_) * (1 - dist_->cdf(x));
  }
};

double weighted_rmst(Distribution * dist, double t, double r, double time_length,
                StateValues * state_values){
  WeightedRmstFunc weighted_rmst_fun(dist, r, time_length, state_values);
  const double lower = 0, upper = t;
  double err_est; int err_code;
  return Numer::integrate(weighted_rmst_fun, lower, upper, err_est, err_code);
}

struct DisModSurv {
  std::string dist_name_;
  vecmats coefs_;
  vecmats X_;
  double time_length_;
  DisModSurv(Rcpp::Environment R_DisModSurv);
  static double summary1(double x, std::string type, Distribution * dist, 
                         double discount_rate, double time_length);
  arma::mat summary(std::vector<double> xvec, std::string type, double discount_rate);
};

class DecModSurv {
private:
  DisModSurv dis_mod_surv_;
  int nsims_; int nindivs_; int npars_;
  
public:
  DecModSurv(DisModSurv dis_mod_surv);
  template <typename Func>
  arma::mat sim_effects(int max_t, Func utility, 
                std::vector<double> discount_rate, std::vector<int> type);
  template <typename Func>
  arma::mat sim_costs(std::vector<int> t, Func cost_values, int n_components,
                std::vector<double> discount_rate);
};

class PartitionSurvival {
private:
  std::string dist_name_;
  vecmats pfs_coefs_;
  vecmats pfs_X_;
  vecmats os_coefs_;
  vecmats os_X_;
  int nsims_; int nindivs_; int npars_;
  double r_health_; double r_costs_; double time_length_;
  
public:
  PartitionSurvival(std::string dist_name, Rcpp::List pfs_coefs, Rcpp::List pfs_X,
                    Rcpp::List os_coefs, Rcpp::List os_X,
                    double r_health, double r_costs, double time_length);
  std::vector<double> predict_pars(int i, int s, int k);
  arma::mat summary(std::vector<double> xvec, std::string type);
};

# endif


