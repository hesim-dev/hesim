#include <hesim/Params.h>

/*************
* Linear model
*************/
ParamsLm::ParamsLm(Rcpp::List R_ParamsLm){
  coefs_ = Rcpp::as<arma::mat>(R_ParamsLm["coefs"]);
  sigma_ = Rcpp::as<std::vector<double> >(R_ParamsLm["sigma"]);
  sample_ = 0;
  n_samples_ = Rcpp::as<int> (R_ParamsLm["n_samples"]);
};

arma::rowvec ParamsLm::get_coefs() const{
  return coefs_.row(sample_);
}

/****************
* Survival models
****************/
SplineSurvAux::SplineSurvAux(Rcpp::List R_ParamsSurv){
  std::string dist_name = Rcpp::as<std::string>(R_ParamsSurv["dist"]);
  if (dist_name == "survspline"){
    Rcpp::List aux = Rcpp::as<Rcpp::List> (R_ParamsSurv["aux"]);
    knots_ = Rcpp::as<std::vector<double> > (aux["knots"]);
    scale_ = Rcpp::as<std::string> (aux["scale"]);
    timescale_  = Rcpp::as<std::string> (aux["timescale"]); 
  }
}

FracPolyAux::FracPolyAux(Rcpp::List R_ParamsSurv){
  std::string dist_name = Rcpp::as<std::string>(R_ParamsSurv["dist"]);
  if (dist_name == "fracpoly"){
    Rcpp::List aux = Rcpp::as<Rcpp::List> (R_ParamsSurv["aux"]);
    powers_ = Rcpp::as<std::vector<double> > (aux["powers"]);
  }
}

// Single survival model
ParamsSurv::ParamsSurv(Rcpp::List R_ParamsSurv)
  : spline_aux_(R_ParamsSurv), fracpoly_aux_(R_ParamsSurv) {
  coefs_ = list_to_vecmats(R_ParamsSurv["coefs"]);
  dist_name_ = Rcpp::as<std::string>(R_ParamsSurv["dist"]);
  sample_ = 0;
  n_samples_ = Rcpp::as<int> (R_ParamsSurv["n_samples"]);
  n_pars_ = coefs_.size();
};

ParamsSurv::ParamsSurv(vecmats coefs, std::string dist_name, int n_samples,
           SplineSurvAux spline_aux, FracPolyAux fracpoly_aux)
  : spline_aux_(spline_aux), fracpoly_aux_(fracpoly_aux){
  coefs_ = coefs;
  dist_name_ = dist_name;
  n_samples_ = n_samples;
};

// List of survival models
ParamsSurvList::ParamsSurvList(Rcpp::List R_ParamsSurvList){
  n_models_ = R_ParamsSurvList.size();
  params_list_.reserve(n_models_);
  for (int i = 0; i < n_models_; ++i){
    Rcpp::List params_surv_i = Rcpp::as<Rcpp::List> (R_ParamsSurvList[i]);
    params_list_.push_back(ParamsSurv(params_surv_i));
  }
  n_samples_ = params_list_[0].n_samples_;
};

// Joined survival models
ParamsJoinedSurvs::ParamsJoinedSurvs(Rcpp::List R_ParamsJoinedSurvs){
  Rcpp::List R_params = Rcpp::as<Rcpp::List> (R_ParamsJoinedSurvs["params"]);
  times_ = Rcpp::as<std::vector<double> > (R_ParamsJoinedSurvs["times"]);
  n_times_ = times_.size();
  for (int i = 0; i < R_ParamsJoinedSurvs.size(); ++i){
    Rcpp::List R_params_i = R_params[i];
    params_[i] = ParamsSurv(R_params_i);
  }
  n_samples_ = params_[0].n_samples_;
}