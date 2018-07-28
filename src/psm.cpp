#include <hesim/psm.h>
#include <hesim/statevals.h>


namespace hesim {

namespace psm {

/****************
* Survival models
****************/
// Base case
surv_mods::surv_mods(Rcpp::Environment R_PartSurvCurves){
  Rcpp::Environment R_data = Rcpp::as<Rcpp::Environment > (R_PartSurvCurves["data"]);
  strategy_id_ = Rcpp::as<std::vector<int> >(R_data["strategy_id"]);
  patient_id_ = Rcpp::as<std::vector<int> >(R_data["patient_id"]);
}

// N separate survival models
surv_list::surv_list(Rcpp::Environment R_PartSurvCurves)
  : surv_mods(R_PartSurvCurves),
    params_(Rcpp::as<Rcpp::List>(R_PartSurvCurves["params"])){
    Rcpp::List R_data = Rcpp::as<Rcpp::List > (R_PartSurvCurves["data"]);
    X_ = Rcpp::as<vecmats_2d>(R_data["X"]);
}

int surv_list::get_n_samples() const {
  return params_.n_samples_;
}

int surv_list::get_n_models() const {
  return params_.n_models_;
}

int surv_list::get_n_obs() const {
  return X_[0][0].n_rows;
}

std::vector<double> surv_list::summary(int model, int sample, int obs, std::vector<double> t,
                                    std::string type, double dr) const{
  hesim::statmods::surv statmod_surv(X_[model], params_.params_list_[model]);
  if (type == "hazard"){
    return statmod_surv.summary(sample, obs, t, "hazard");
  }
  else if (type == "cumhazard"){
    return statmod_surv.summary(sample, obs, t, "cumhazard");
  }
  else if (type == "survival"){
    return statmod_surv.summary(sample, obs, t, "survival");
  }
  else if (type == "rmst"){
    return statmod_surv.summary(sample, obs, t, "rmst", dr);
  }
  else {
    Rcpp::stop("Selected type is not available.");
  }
}

std::vector<double> surv_list::quantile(int model, int sample, int obs, std::vector<double> p) const{
  hesim::statmods::surv statmod_surv(X_[model], params_.params_list_[model]);
  return statmod_surv.quantile(sample, obs, p);
}

// Create pointer to base survival class
std::unique_ptr<surv_mods> create_surv_mods(Rcpp::Environment R_PartSurvCurves){
    Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_PartSurvCurves["params"]);
 
    surv_mods * survmods;
    if (Rf_inherits(R_params, "params_surv_list")){
        survmods = new surv_list(R_PartSurvCurves);
    }
    else{
        Rcpp::stop("The selected statistical model is not available.");
    }
    std::unique_ptr<surv_mods> uptr(survmods);
    return uptr;
}

/****************
* Survival curves
****************/
surv_curves_out::surv_curves_out(){
}

surv_curves_out::surv_curves_out(int n){
  curve_.resize(n);
  sample_.resize(n);
  strategy_id_.resize(n);
  patient_id_.resize(n);
  x_.resize(n);
  value_.resize(n);
}

int surv_curves_out::index(int sim, int curve) const{
  return curve * n_sims_ + sim;
}

surv_curves_out surv_curves_out::create_from_R(Rcpp::Environment R_PartSurv){
  surv_curves_out out;
  Rcpp::DataFrame R_survival = Rcpp::as<Rcpp::DataFrame>(R_PartSurv["survival_"]);
  
  out.curve_ = Rcpp::as<std::vector<int> >(R_survival["curve"]);
  out.sample_ = Rcpp::as<std::vector<int> >(R_survival["sample"]);
  out.strategy_id_ = Rcpp::as<std::vector<int> >(R_survival["strategy_id"]);
  out.patient_id_ = Rcpp::as<std::vector<int> >(R_survival["patient_id"]);
  out.x_ = Rcpp::as<std::vector<double> >(R_survival["t"]);
  out.value_ = Rcpp::as<std::vector<double> >(R_survival["survival"]);
  out.n_states_ = Rcpp::as<int>(R_PartSurv["n_states"]);
  int n_curves = out.n_states_ - 1;
  out.n_sims_ = out.value_.size()/n_curves;
  
  // R to C++ indexing
  hesim::add_constant(out.curve_, -1); 
  hesim::add_constant(out.sample_, -1); 
  return out;
}

surv_curves::surv_curves(Rcpp::Environment R_PartSurvCurves)
  : survmods_(create_surv_mods(R_PartSurvCurves)){
};

Rcpp::DataFrame surv_curves::summary(std::vector<double> x, std::string type,
                                double dr){
  // Preallocate
  int n_models = survmods_->get_n_models();
  int n_samples = survmods_->get_n_samples();
  int n_obs = survmods_->get_n_obs();
  int N = n_models * n_samples * n_obs * x.size();
  surv_curves_out out(N);            

  // Loop
  int counter = 0;
  for (int m = 0; m < n_models; ++m){
    for (int s = 0; s < n_samples; ++s){
      for (int i = 0; i < n_obs; ++i){
        std::vector<double> res_vec;
        if (type == "survival"){
          res_vec = survmods_->summary(m, s, i, x, "survival");
        }
        else if (type == "cumhazard"){
          res_vec = survmods_->summary(m, s, i, x, "cumhazard");
        }
        else if (type == "hazard"){
          res_vec = survmods_->summary(m, s, i, x, "hazard");
        }
        else if (type == "rmst"){
          res_vec = survmods_->summary(m, s, i, x, "rmst", dr); 
        }
        else if (type == "quantile"){
          res_vec = survmods_->quantile(m, s, i, x);
        } 
        else{
          Rcpp::stop("Selected type for summarizing survival distribution is not available.");
        }
        
        // store vector
        for (int j = 0; j < x.size(); ++j){
          out.curve_[counter] = m;
          out.sample_[counter] = s;
          out.strategy_id_[counter] = survmods_->strategy_id_[i];
          out.patient_id_[counter] = survmods_->patient_id_[i];
          out.x_[counter] = x[j];
          out.value_[counter] = res_vec[j];
          ++counter;
        }
      }
    }
  }

  // Return
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["curve"] = out.curve_,
    Rcpp::_["sample"] = out.sample_,
    Rcpp::_["strategy_id"] = out.strategy_id_,
    Rcpp::_["patient_id"] = out.patient_id_,
    Rcpp::_["x"] = out.x_,
    Rcpp::_["value"] = out.value_,
    Rcpp::_["stringsAsFactors"] = false
  );
  return out_df;
}

/************************************
* Partitioned survival decision model
************************************/
stateprobs::stateprobs(Rcpp::Environment R_PartSurv)
  : survcurves_(surv_curves_out::create_from_R(R_PartSurv)){
  n_crossings_ = 0; // Number of times survival curves cross
};

double stateprobs::sim_probs1(int sim, int state){
  if (state == 0){
    int index = survcurves_.index(sim, state);
    return survcurves_.value_[index];
  }
  else if (state > 0 && state < (survcurves_.n_states_ - 1)){
    int index = survcurves_.index(sim, state);
    int index_prior = survcurves_.index(sim, state - 1);
    double survival = survcurves_.value_[index];
    double survival_prior = survcurves_.value_[index_prior];
    if (survival < survival_prior){
      survival = survival_prior;
      ++n_crossings_;
    }
    return survival - survival_prior;
  }
  else {
    int index_prior = survcurves_.index(sim, state - 1);
    return 1 - survcurves_.value_[index_prior];
  }
}

Rcpp::List stateprobs::sim_probs(){
  stateprobs_out out(survcurves_.n_sims_ * survcurves_.n_states_ );

  int counter = 0;
  for (int j = 0; j < survcurves_.n_states_; ++j){
    for (int i = 0; i < survcurves_.n_sims_; ++i){
      int index = survcurves_.index(i, 0);
      out.state_id_[counter] = j;
      out.sample_[counter] = survcurves_.sample_[index];
      out.strategy_id_[counter] = survcurves_.strategy_id_[index];
      out.patient_id_[counter] = survcurves_.patient_id_[index];
      out.t_[counter] = survcurves_.x_[index];
      out.prob_[counter] = sim_probs1(i, j);
      ++counter;
    }
  }
  
  // Return
  Rcpp::DataFrame stateprobs_df = Rcpp::DataFrame::create(
    Rcpp::_["state_id"] = out.state_id_,
    Rcpp::_["sample"] = out.sample_,
    Rcpp::_["strategy_id"] = out.strategy_id_,
    Rcpp::_["patient_id"] = out.patient_id_,
    Rcpp::_["t"] = out.t_,
    Rcpp::_["prob"] = out.prob_,
    Rcpp::_["stringsAsFactors"] = false
  );
  
  return(Rcpp::List::create(
    Rcpp::_["stateprobs"] = stateprobs_df,
    Rcpp::_["n_crossings"] = n_crossings_
  ));
};

} // end psm namespace

} // end hesim namespace

/************************
* Functions exported to R
************************/
// [[Rcpp::export]]
Rcpp::DataFrame C_PartSurvCurves_summary(Rcpp::Environment R_PartSurvCurves, std::vector<double> x, 
                             std::string type, double dr){
  hesim::psm::surv_curves survcurves(R_PartSurvCurves);
  return survcurves.summary(x, type, dr);
}

// [[Rcpp::export]]
Rcpp::List C_PartSurv_sim_stateprobs(Rcpp::Environment R_PartSurv){
  hesim::psm::stateprobs stprobs(R_PartSurv);
  return stprobs.sim_probs();
}

// [[Rcpp::export]]
Rcpp::DataFrame C_psm_sim_wlos(Rcpp::Environment R_psm, Rcpp::DataFrame R_stateprobs, 
                              std::vector<double> dr, std::string type,
                              std::vector<std::string> categories){
  hesim::wlos wlos(R_psm, type);
  hesim::stateprobs_out stprobs(R_stateprobs);
  std::vector<double> times = Rcpp::as<std::vector<double> > (R_psm["t_"]);
  hesim::wlos_out out = wlos(stprobs, times, dr, categories);
  return out.create_R_data_frame();
}




