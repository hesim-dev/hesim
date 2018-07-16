#include "PartSurv.h"

namespace part_surv{

/*************************
* Models for state values
*************************/
// Base case
StateValMods::StateValMods(Rcpp::Environment R_PartSurvStateVals){
  Rcpp::Environment R_data = Rcpp::as<Rcpp::Environment > (R_PartSurvStateVals["data"]);
  state_id_ = Rcpp::as<std::vector<int> >(R_data["state_id"]);
  strategy_id_ = Rcpp::as<std::vector<int> >(R_data["strategy_id"]);
  patient_id_ = Rcpp::as<std::vector<int> >(R_data["patient_id"]);
  n_strategies_ = Rcpp::as<int>(R_data["n_strategies"]);
  n_patients_ = Rcpp::as<int>(R_data["n_patients"]);
  n_states_ = Rcpp::as<int>(R_data["n_states"]);
  obs_ = 0;
  sample_ = 0;
}

void StateValMods::set_sample(int sample){
  sample_ = sample;
}

// Linear model
LmMod::LmMod(Rcpp::Environment R_PartSurvStateVals)
  : StateValMods(R_PartSurvStateVals),
    params_(Rcpp::as<Rcpp::List>(R_PartSurvStateVals["params"])){
  Rcpp::List R_data = Rcpp::as<Rcpp::List > (R_PartSurvStateVals["data"]);
  Rcpp::List R_X = Rcpp::as<Rcpp::List > (R_data["X"]);
  X_ = Rcpp::as<arma::mat>(R_X["mu"]);
}

void LmMod::set_obs(int strategy, int patient, int state){
  obs_ = strategy * n_patients_ * n_states_ +
         patient * n_states_ +
         state;
} 

void LmMod::set_obs(int obs){
  obs_ = obs;
}

int LmMod::get_n_obs(){
  return X_.n_rows;
}

int LmMod::get_n_samples(){
  return params_.n_samples_;
}

double LmMod::predict() const{
  StatModLm statmod_lm(X_, params_);
  return statmod_lm.predict(sample_, obs_);
}

// Create pointer to base state values model class
std::unique_ptr<StateValMods> create_StateValMods(Rcpp::Environment R_PartSurvStateVal){
    Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_PartSurvStateVal["params"]);
 
    StateValMods * statevalmods;
    if (Rf_inherits(R_params, "params_lm")){
        statevalmods = new LmMod(R_PartSurvStateVal);
    }
    else{
        Rcpp::stop("The selected statistical model is not available.");
    }
    std::unique_ptr<StateValMods> uptr(statevalmods);
    return uptr;
}


/****************
* Survival models
****************/
// Base case
SurvMods::SurvMods(Rcpp::Environment R_PartSurvCurves){
  Rcpp::Environment R_data = Rcpp::as<Rcpp::Environment > (R_PartSurvCurves["data"]);
  strategy_id_ = Rcpp::as<std::vector<int> >(R_data["strategy_id"]);
  patient_id_ = Rcpp::as<std::vector<int> >(R_data["patient_id"]);
}

// N separate survival models
NSurvMods::NSurvMods(Rcpp::Environment R_PartSurvCurves)
  : SurvMods(R_PartSurvCurves),
    params_(Rcpp::as<Rcpp::List>(R_PartSurvCurves["params"])){
    Rcpp::List R_data = Rcpp::as<Rcpp::List > (R_PartSurvCurves["data"]);
    X_ = Rcpp::as<vecmats_2d>(R_data["X"]);
}

int NSurvMods::get_n_samples() const {
  return params_.n_samples_;
}

int NSurvMods::get_n_models() const {
  return params_.n_models_;
}

int NSurvMods::get_n_obs() const {
  return X_[0][0].n_rows;
}

std::vector<double> NSurvMods::summary(int model, int sample, int obs, std::vector<double> t,
                                    std::string type, double dr) const{
  StatModSurv statmod_surv(X_[model], params_.params_list_[model]);
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

std::vector<double> NSurvMods::quantile(int model, int sample, int obs, std::vector<double> p) const{
  StatModSurv statmod_surv(X_[model], params_.params_list_[model]);
  return statmod_surv.quantile(sample, obs, p);
}

// // N separate joined survival models
// NJoinedSurvs::NJoinedSurvs(Rcpp::Environment R_Curves)
//   : Survs(R_Curves),
//     params_(Rcpp::as<Rcpp::Environment>(R_Curves["params"])){
//   
//   Rcpp::Environment R_data = Rcpp::as<Rcpp::Environment > (R_Curves["data"]);
//   X_ = Rcpp::as<vecmats_3d>(R_data["X"]);
// }
// 
// std::vector<double> NJoinedSurvs::predict_pars() const{
//   int n_pars = params_.coefs_[model_][time_id_].size();
//   std::vector<double> y(n_pars);
//   for (int j = 0; j < n_pars; ++j){
//     y[j] = arma::dot(X_[model_][time_id_][j].row(obs_), 
//                      params_.coefs_[model_][time_id_][j].row(sim_));
//   }
//   return y;
// }
// 
// void NJoinedSurvs::set_dist(){
//   dist_ = hesim::select_distribution(params_.dist_names_[model_][time_id_],
//                                      predict_pars());
// }
// 
// int NJoinedSurvs::get_n_sims() const {
//   return params_.n_sims_;
// }
// 
// int NJoinedSurvs::get_n_models() const {
//   return params_.n_models_;
// }
// 
// void NJoinedSurvs::time_update(double t){
//   if (t > params_.times[model_][time_id_]){
//     ++time_id_;
//     set_dist();
//   }
// }
// 
// template <typename T>
// std::vector<double> NJoinedSurvs::outcomes(std::vector<double> x, T fun){
//   std::vector<double> out;
//   time_id_ = 0; 
//   for (int i = 0; i < x.size(); ++i){
//     time_update(x[i]);
//     out[i] = fun(x[i]);
//   }
//   return out;
// }
// 
// std::vector<double> NJoinedSurvs::hazard(std::vector<double> t){
//   hesim::HazardFunc hazard(dist_);
//   return outcomes(t, hazard);
// }
// 
// std::vector<double> NJoinedSurvs::cumhazard(std::vector<double> t){
//   hesim::CumhazardNumericFunc cumhazard(dist_);
//   return outcomes(t, cumhazard);
// }
// 
// std::vector<double> NJoinedSurvs::survival(std::vector<double> t){
//   hesim::SurvivalNumericFunc surv(dist_);
//   return outcomes(t, surv);
// }
// 
// std::vector<double> NJoinedSurvs::quantiles(std::vector<double> p){
//   hesim::QuantileNumericFunc quants(dist_);
//   return outcomes(p, quants);
// }

// Create pointer to base survival class
std::unique_ptr<SurvMods> create_SurvMods(Rcpp::Environment R_PartSurvCurves){
    Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_PartSurvCurves["params"]);
 
    SurvMods * survmods;
    if (Rf_inherits(R_params, "params_surv_list")){
        survmods = new NSurvMods(R_PartSurvCurves);
    }
    else{
        Rcpp::stop("The selected statistical model is not available.");
    }
    std::unique_ptr<SurvMods> uptr(survmods);
    return uptr;
}

/****************
* Survival curves
****************/
SurvCurvesOut::SurvCurvesOut(){
}

SurvCurvesOut::SurvCurvesOut(int n){
  curve_.resize(n);
  sample_.resize(n);
  strategy_id_.resize(n);
  patient_id_.resize(n);
  x_.resize(n);
  value_.resize(n);
}

int SurvCurvesOut::index(int sim, int curve) const{
  return curve * n_sims_ + sim;
}

SurvCurvesOut SurvCurvesOut::create_from_R(Rcpp::Environment R_PartSurv){
  SurvCurvesOut out;
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

SurvCurves::SurvCurves(Rcpp::Environment R_PartSurvCurves)
  : survmods_(create_SurvMods(R_PartSurvCurves)){
};

Rcpp::DataFrame SurvCurves::summary(std::vector<double> x, std::string type,
                                double dr){
  // Preallocate
  int n_models = survmods_->get_n_models();
  int n_samples = survmods_->get_n_samples();
  int n_obs = survmods_->get_n_obs();
  int N = n_models * n_samples * n_obs * x.size();
  SurvCurvesOut out(N);            

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

/*************
* State values
*************/
StateValsOut::StateValsOut(int n){
  state_id_.resize(n);
  sample_.resize(n);
  strategy_id_.resize(n);
  patient_id_.resize(n);
  value_.resize(n);
}

StateVals::StateVals(Rcpp::Environment R_PartSurvStateVals)
  : statevalmods_(create_StateValMods(R_PartSurvStateVals)){
};

Rcpp::DataFrame StateVals::predict(){
  int n_obs = statevalmods_->get_n_obs();
  int n_samples = statevalmods_->get_n_samples();
  StateValsOut out(n_obs * n_samples);

  int counter = 0;
  for (int s = 0; s < n_samples; ++ s){
    statevalmods_->set_sample(s);
    for (int i = 0; i < n_obs; ++i){
      out.state_id_[counter] = statevalmods_->state_id_[i];
      out.sample_[counter] = statevalmods_->sample_;
      out.strategy_id_[counter] = statevalmods_->strategy_id_[i];
      out.patient_id_[counter] = statevalmods_->patient_id_[i];
      statevalmods_->set_obs(i);
      out.value_[counter] = statevalmods_->predict();
      ++ counter;
    }
  }
  
  // Return
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["state_id"] = out.state_id_,
    Rcpp::_["sample"] = out.sample_,
    Rcpp::_["strategy_id"] = out.strategy_id_,
    Rcpp::_["patient_id"] = out.patient_id_,
    Rcpp::_["value"] = out.value_,
    Rcpp::_["stringsAsFactors"] = false
  );
  return out_df;
};

/************************************
* Partitioned survival decision model
************************************/
StateprobsOut::StateprobsOut(){
}

StateprobsOut::StateprobsOut(int n){
  state_id_.resize(n);
  sample_.resize(n);
  strategy_id_.resize(n);
  patient_id_.resize(n);
  t_.resize(n);
  prob_.resize(n);
}

StateprobsOut StateprobsOut::create_from_R(Rcpp::DataFrame R_stateprobs){
  StateprobsOut out;
  
  out.state_id_ = Rcpp::as<std::vector<int> >(R_stateprobs["state_id"]);
  out.sample_ = Rcpp::as<std::vector<int> >(R_stateprobs["sample"]);
  out.strategy_id_ = Rcpp::as<std::vector<int> >(R_stateprobs["strategy_id"]);
  out.patient_id_ = Rcpp::as<std::vector<int> >(R_stateprobs["patient_id"]);
  out.t_ = Rcpp::as<std::vector<double> >(R_stateprobs["t"]);
  out.prob_ = Rcpp::as<std::vector<double> >(R_stateprobs["prob"]);

  // R to C++ indexing
  hesim::add_constant(out.state_id_, -1);
  hesim::add_constant(out.sample_, -1);
  hesim::add_constant(out.strategy_id_, -1);
  hesim::add_constant(out.patient_id_, -1);
  return out;
}

Stateprobs::Stateprobs(Rcpp::Environment R_PartSurv)
  : survcurves_(SurvCurvesOut::create_from_R(R_PartSurv)){
  n_crossings_ = 0; // Number of times survival curves cross
};

double Stateprobs::sim_probs1(int sim, int state){
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

Rcpp::List Stateprobs::sim_probs(){
  StateprobsOut out(survcurves_.n_sims_ * survcurves_.n_states_ );

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

AucOut::AucOut(int n){
  state_id_.resize(n);
  sample_.resize(n);
  strategy_id_.resize(n);
  patient_id_.resize(n);
  dr_.resize(n);
  type_.resize(n);
  value_.resize(n);
}

Auc::Auc(Rcpp::Environment R_PartSurv, Rcpp::DataFrame R_stateprobs, std::string type)
  : stateprobs_(StateprobsOut::create_from_R(R_stateprobs)),
    statevalmods_(create_statevalmods_(R_PartSurv, type)){
  times_ = Rcpp::as<std::vector<double> > (R_PartSurv["t_"]);
  n_times_ = times_.size();
  n_states_ = Rcpp::as<int>(R_PartSurv["n_states"]);
}

std::vector<std::unique_ptr<StateValMods> > Auc::create_statevalmods_(Rcpp::Environment R_PartSurv, 
                                                        std::string type){
  std::vector<std::unique_ptr<StateValMods> > statevalmods;
  if (type == "qalys"){
    Rcpp::Environment R_utility_model = Rcpp::as<Rcpp::Environment>(R_PartSurv["utility_model"]);
    statevalmods.push_back(create_StateValMods(R_utility_model));
  } 
  else if (type == "costs"){
    Rcpp::List R_costs_models = Rcpp::as<Rcpp::List>((R_PartSurv["cost_models"]));
    for (int i = 0; i < R_costs_models.size(); ++i){
      Rcpp::Environment R_cost_model_i = Rcpp::as<Rcpp::Environment>(R_costs_models[i]);
      statevalmods.push_back(create_StateValMods(R_cost_model_i));
    }
  }
  else{
    Rcpp::stop("Predictions can only be made for 'costs' or 'qalys'.");
  }
  return statevalmods;
}

// Calculated weighted area under survival curve using trapezoid method
double auc_trapz(std::vector<double>::iterator t_first, std::vector<double>::iterator t_last,
                   std::vector<double>::iterator stateprob_first, StateValMods *statevalmods, 
                   double dr){
  std::vector<double> value(std::distance(t_first, t_last));
  std::vector<double>::iterator stateprob_it = stateprob_first;
  for (std::vector<double>::iterator t_it = t_first; t_it != t_last; ++t_it){
    double statval = statevalmods->predict();
    value[t_it - t_first] = exp(-dr * *t_it) * statval * *stateprob_it;
    ++stateprob_it;
  }
  return trapz(t_first, t_last, value.begin());
}

Rcpp::DataFrame Auc::sim(std::vector<double> dr, std::vector<std::string> type_names){
  int N = stateprobs_.prob_.size()/n_times_;
  AucOut out(N * dr.size() * statevalmods_.size());
  
  int counter = 0;
  for (int k = 0; k < statevalmods_.size(); ++k){
    for (int j = 0; j < dr.size(); ++j){
      int index = 0;
      double dr_j = dr[j];
      for (int i = 0; i < N; ++i){
        int state_id = stateprobs_.state_id_[index];
        int sample = stateprobs_.sample_[index];
        int strategy_id = stateprobs_.strategy_id_[index];
        int patient_id = stateprobs_.patient_id_[index];
        
        out.state_id_[counter] = state_id;
        out.sample_[counter] = sample;
        out.strategy_id_[counter] = strategy_id;
        out.patient_id_[counter] = patient_id;
        out.dr_[counter] = dr_j;
        out.type_[counter] = type_names[k];
        statevalmods_[k]->set_sample(stateprobs_.sample_[index]);
        statevalmods_[k]->set_obs(strategy_id, patient_id, state_id);
        out.value_[counter] = auc_trapz(times_.begin(), times_.end(),
                                         stateprobs_.prob_.begin() + index, 
                                         statevalmods_[k].get(), 
                                         dr_j);
        index = index + times_.size();
        ++counter;
      } // end loop over state probabilities
    } // end loop over discount rates
  } // end loop over state value models

  // Return
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["state_id"] = out.state_id_,
    Rcpp::_["sample"] = out.sample_,
    Rcpp::_["strategy_id"] = out.strategy_id_,
    Rcpp::_["patient_id"] = out.patient_id_,
    Rcpp::_["dr"] = out.dr_,
    Rcpp::_["type"] = out.type_,
    Rcpp::_["value"] = out.value_,
    Rcpp::_["stringsAsFactors"] = false
  );
  
  return out_df;
}

} // end part_surv namespace

/************************
* Functions exported to R
************************/
// [[Rcpp::export]]
Rcpp::DataFrame C_PartSurvCurves_summary(Rcpp::Environment R_PartSurvCurves, std::vector<double> x, 
                             std::string type, double dr){
  part_surv::SurvCurves survcurves(R_PartSurvCurves);
  return survcurves.summary(x, type, dr);
}

// [[Rcpp::export]]
Rcpp::DataFrame C_PartSurvStateVals_predict(Rcpp::Environment R_PartSurvStateVals){
  part_surv::StateVals statevals(R_PartSurvStateVals);
  return statevals.predict();
}

// [[Rcpp::export]]
Rcpp::List C_PartSurv_sim_stateprobs(Rcpp::Environment R_PartSurv){
  part_surv::Stateprobs stateprobs(R_PartSurv);
  return stateprobs.sim_probs();
}

// [[Rcpp::export]]
Rcpp::List C_PartSurv_sim_auc(Rcpp::Environment R_PartSurv, Rcpp::DataFrame R_stateprobs, 
                              std::vector<double> dr, std::string type,
                              std::vector<std::string> type_names){
  part_surv::Auc auc(R_PartSurv, R_stateprobs, type);
  return auc.sim(dr, type_names);
}


