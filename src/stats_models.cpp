#include <hesim/stats_models.h>

/*******************
* Linear model (OLS)
*******************/
StatModLm::StatModLm(arma::mat X, ParamsLm params)
  : params_(params){
  X_ = X;
}

double StatModLm::predict(int sample, int obs){
  return arma::dot(X_.row(obs), params_.coefs_.row(sample));
}

double StatModLm::simulate(int sample, int obs){
  return R::rnorm(predict(sample, obs), params_.sigma_[sample]);
}

/****************
* Survival models
****************/
// Initialize distribution
std::unique_ptr<hesim::Distribution> init_Distribution_ptr(ParamsSurv params_surv){
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

// Survival model
std::vector<double> StatModSurv::predict_pars(int sample, int obs) const{
  int n_pars = params_.coefs_.size();
  std::vector<double> y(n_pars);
  for (int j = 0; j < n_pars; ++j){
    y[j] = arma::dot(X_[j].row(obs),
                     params_.coefs_[j].row(sample));
  }
  return y;
};

StatModSurv::StatModSurv(vecmats X, ParamsSurv params)
  : params_(params), dist_(init_Distribution_ptr(params)){
  X_ = X;
}

void StatModSurv::set_dist(int sample, int obs){
  dist_->set_params(predict_pars(sample, obs));
}

std::vector<double> StatModSurv::summary(int sample, int obs, std::vector<double> t, std::string type,
                            double dr){
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

std::vector<double> StatModSurv::quantile(int sample, int obs, std::vector<double> p){
  std::vector<double> out(p.size());
  set_dist(sample, obs);
  for (int i = 0; i < p.size(); ++i){
    out[i] = dist_->quantile(p[i]);
  }
  return out;
};