# ifndef STATMODELS_H
# define STATMODELS_H
#include <hesim/distributions.h>
#include "Params.h"

/*******************
* Linear model (OLS)
*******************/
class StatModLm {
public:
  arma::mat X_;
  ParamsLm params_;
  StatModLm(arma::mat X, ParamsLm params)
    : params_(params){
    X_ = X;
  }
  
  double predict(int sample, int obs){
    return arma::dot(X_.row(obs), params_.coefs_.row(sample));
  }
  
  double simulate(int sample, int obs){
    return R::rnorm(predict(sample, obs), params_.sigma_[sample]);
  }
};

/****************
* Survival models
****************/
// Initialize distribution pointer with dummy values for parameters and actual values for
// auxillary parameters
inline std::unique_ptr<hesim::Distribution> init_Distribution_ptr(ParamsSurv params_surv){
  hesim::Distribution *d;
  std::string dist_name = params_surv.dist_name_;
  if (dist_name == "exponential" || dist_name == "exp"){
    d = new hesim::Exponential(1);
  }
  else if (dist_name == "weibull.quiet" || dist_name == "weibull"){
    d = new hesim::Weibull(1, 1);
  }
  else if (dist_name == "weibullNMA"){
    d = new hesim::WeibullNma(0, 0); // equivalent to shape = 1, scale = 1 wih weibull
  }
  else if (dist_name == "gamma"){
    d = new hesim::Gamma(1, 1);
  }
  else if (dist_name == "lnorm"){
    d = new hesim::Lognormal(0, 1);
  }
  else if (dist_name == "gompertz"){
    d = new hesim::Gompertz(0, 1);
  }
  else if (dist_name == "llogis"){
    d = new hesim::LogLogistic(1, 1);
  }
  else if (dist_name == "gengamma"){
    d = new hesim::GeneralizedGamma(0, 1, 0);
  }
  else if (dist_name == "survspline"){
    int n_knots = params_surv.spline_aux_.knots_.size();
    std::vector<double> gamma = zeros<double>(n_knots);
    d = new hesim::SurvSplines(gamma, params_surv.spline_aux_.knots_,
                               params_surv.spline_aux_.scale_,
                               params_surv.spline_aux_.timescale_);
  }
  else{
      Rcpp::stop("The selected distribution is not available.");
  }
  std::unique_ptr<hesim::Distribution> uptr(d);
  return uptr;
}


class StatModSurv {
private:
  std::vector<double> predict_pars(int sample, int obs) const{
    int n_pars = params_.coefs_.size();
    std::vector<double> y(n_pars);
    for (int j = 0; j < n_pars; ++j){
      y[j] = arma::dot(X_[j].row(obs),
                       params_.coefs_[j].row(sample));
    }
    return y;
  };

public:
  vecmats X_;
  ParamsSurv params_;
  std::unique_ptr<hesim::Distribution> dist_;
  StatModSurv(vecmats X, ParamsSurv params)
    : params_(params), dist_(init_Distribution_ptr(params)){
    X_ = X;
  }

  void set_dist(int sample, int obs){
    dist_->set_params(predict_pars(sample, obs));
  }

  // Return fitted survival, cummulative hazard, hazard, or restricted mean
  // survival time at a series of times
  std::vector<double> summary(int sample, int obs, std::vector<double> t, std::string type,
                              double dr = 0){
    std::vector<double> out(t.size());
    set_dist(sample, obs);
    for (int i = 0; i < t.size(); ++i){
      if (type == "hazard"){
        out[i] = dist_->hazard(t[i]);
      }
      else if (type == "cumhazard"){
        out[i] = dist_->cumhazard(t[i]);
      }
      else if (type == "survival"){
        out[i] = 1 - dist_->cdf(t[i]);
      }
      else if (type == "rmst"){
        out[i] = hesim::rmst(dist_.get(), t[i], dr);
      }
      else{
        Rcpp::stop("Selected 'type' is not supported.");
      }
    }
    return out;
  };

  // Return quantile of distribution
  std::vector<double> quantile(int sample, int obs, std::vector<double> p){
    std::vector<double> out(p.size());
    set_dist(sample, obs);
    for (int i = 0; i < p.size(); ++i){
      out[i] = dist_->quantile(p[i]);
    }
    return out;
  };

};

// Initialize a vector of StatModSurv objects given a sampled parameter set
inline std::vector<StatModSurv> init_joined_statmod_survs(int sample, int obs,
                                                          vecmats_2d X,
                                                          ParamsJoinedSurvs params){
  int N = params.params_.size();
  std::vector<StatModSurv> joined_statmod_survs;
  for (int i = 0; i < N; ++i){
    StatModSurv statmod_surv_i(X[i], params.params_[i]);
    statmod_surv_i.set_dist(sample, obs);
    joined_statmod_survs.push_back(std::move(statmod_surv_i)); //statmod_surv_i contains a unique_ptr
  }
  return joined_statmod_survs;
}

// A class storing the joined survival curves with a member function
// to select the survival model used at time t
class JoinedSurvs {
private:
  int n_times_;
  std::vector<double> times_;
public:
  std::vector<StatModSurv> statmods_;

  JoinedSurvs(int sample, int obs, vecmats_2d X, ParamsJoinedSurvs params)
    : statmods_(init_joined_statmod_survs(sample, obs, X, params)) {
    n_times_ = params.n_times_;
    times_ = params.times_;
  };

  int get_time_id(double t) const{
    for (int i = 0; i < n_times_; ++i){
      if (t < times_[i]){
        return i;
      }
    }
    return n_times_;
  };

  // TO DO: THIS IS MORE EFFICIENT BUT DOESN'T WORK WITH CONST OPERATOR IN RCPP:NUMERICAL
  // PERHAPS REMOVE DEPENDENCY ON RCPP NUMERICAL AND KEEP TIME_ID_ AS A MEMBER VARIABLE.
  //  void update_time (double t) const{
  //   if (t > params_.times_[time_id_]){
  //     return ++time_id_;
  //   }
  // };

};


// Summarize joined survival curves by integrating over hazard function
class JoinedSurvsHazard: public Numer::Func {
public:
  JoinedSurvs joined_survs_;
  JoinedSurvsHazard(int sample, int obs, vecmats_2d X,
                    ParamsJoinedSurvs params)
    : joined_survs_(sample, obs, X, params) {
  };
  double operator()(const double& x) const {
    int time_id = joined_survs_.get_time_id(x);
    return joined_survs_.statmods_[time_id].dist_->hazard(x);
  }
};

class JoinedSurvsRmst: public Numer::Func {
public:
  JoinedSurvs joined_survs_;
  double dr_;
  JoinedSurvsRmst(int sample, int obs, vecmats_2d X,
                    ParamsJoinedSurvs params, double dr)
    : joined_survs_(sample, obs, X, params) {
    dr_ = dr;
  };
  double operator()(const double& x) const {
    int time_id = joined_survs_.get_time_id(x);
    return exp(-dr_ * x) * (1 - joined_survs_.statmods_[time_id].dist_->cdf(x));
  }
};


class StatModJoinedSurvs {
public:
  JoinedSurvs joined_survs_;
  vecmats_2d X_;
  ParamsJoinedSurvs params_;

  StatModJoinedSurvs(vecmats_2d X, ParamsJoinedSurvs params)
    : joined_survs_(0, 0, X, params), params_(params){
    X_ = X;
  }

  double hazard(int sample, int obs, double t){
    JoinedSurvsHazard hazard(sample, obs, X_, params_);
    return hazard(t);
  };

  double cumhazard(int sample, int obs, double t){
    JoinedSurvsHazard hazard(sample, obs, X_, params_);
    const double lower = 0, upper = t;
    double err_est; int err_code;
    return Numer::integrate(hazard, lower, upper, err_est, err_code);
  };

  double survival(int sample, int obs, double t){
    return exp(-cumhazard(sample, obs, t));
  };

  double rmst(int sample, int obs, double t, double dr){
    JoinedSurvsRmst fun(sample, obs, X_, params_, dr);
    const double lower = 0, upper = t;
    double err_est; int err_code;
    return Numer::integrate(fun, lower, upper, err_est, err_code);
  };


  // Return fitted survival, cummulative hazard, hazard, or restricted mean
  // survival time at a series of times
  std::vector<double> summary(int sample, int obs, std::vector<double> t, std::string type,
                              double dr = 0){
    std::vector<double> out(t.size());
    for (int i = 0; i < t.size(); ++i){
      if (type == "hazard"){
        out[i] = hazard(sample, obs, t[i]);
      }
      else if (type == "cumhazard"){
        out[i] = cumhazard(sample, obs, t[i]);
      }
      else if (type == "survival"){
        out[i] = survival(sample, obs, t[i]);
      }
      else if (type == "rmst"){
        out[i] = rmst(sample, obs, t[i], dr);
      }
      else{
        Rcpp::stop("Selected 'type' is not supported.");
      }
    }
    return out;
  };

  // Return quantile of distribution

}; // End StatModJoinedSurvs


# endif


