#include <hesim/psm.h>
#include <hesim/statevals.h>


namespace hesim {

namespace psm {

/****************
* Survival models
****************/
// Base case
surv_mods::surv_mods(Rcpp::Environment R_PsmCurves){
  Rcpp::Environment R_data = Rcpp::as<Rcpp::Environment > (R_PsmCurves["data"]);
  strategy_id_ = Rcpp::as<std::vector<int> >(R_data["strategy_id"]);
  patient_id_ = Rcpp::as<std::vector<int> >(R_data["patient_id"]);
}

std::unique_ptr<surv_mods> surv_mods::create(Rcpp::Environment R_PsmCurves){
    Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_PsmCurves["params"]);
 
    surv_mods * survmods;
    if (Rf_inherits(R_params, "params_surv_list")){
        survmods = new surv_list(R_PsmCurves);
    }
    else{
        Rcpp::stop("The selected statistical model is not available.");
    }
    std::unique_ptr<surv_mods> uptr(survmods);
    return uptr;
}


// N separate survival models
surv_list::surv_list(Rcpp::Environment R_PsmCurves)
  : surv_mods(R_PsmCurves),
    params_(Rcpp::as<Rcpp::List>(R_PsmCurves["params"])){
    Rcpp::List R_data = Rcpp::as<Rcpp::List > (R_PsmCurves["data"]);
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
std::unique_ptr<surv_mods> create_surv_mods(Rcpp::Environment R_PsmCurves){
    Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_PsmCurves["params"]);
 
    surv_mods * survmods;
    if (Rf_inherits(R_params, "params_surv_list")){
        survmods = new surv_list(R_PsmCurves);
    }
    else{
        Rcpp::stop("The selected statistical model is not available.");
    }
    std::unique_ptr<surv_mods> uptr(survmods);
    return uptr;
}

/****************************
* Health state probabilities
****************************/
stateprobs::stateprobs(Rcpp::Environment R_psm)
  : survcurves_(surv_summary::create_from_R(R_psm)){
  n_crossings_ = 0; // Number of times survival curves cross
};

double stateprobs::sim1(int sim, int state){
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

Rcpp::List stateprobs::sim(){
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
      out.prob_[counter] = sim1(i, j);
      ++counter;
    } // end health state loop
  } // end simulation loop
  
  // Return
  Rcpp::DataFrame stateprobs_df = out.create_R_data_frame();
  
  return(Rcpp::List::create(
    Rcpp::_["stateprobs"] = stateprobs_df,
    Rcpp::_["n_crossings"] = n_crossings_
  ));
};

} // end psm namespace

} // end hesim namespace

/***************************************************************************//** 
 * @ingroup psm
 * Summarize the survival curves from a partitioned survival model.
 * This function is exported to @c R and used in the private member function
 *  @c PsmCurves$summary().
 * @param R_psm_curves An R object of class @c PsmCurves.
 * @param x A vector of values at which to summarize the survival curves 
 * (hazard, cumulative hazard, survival, restricted mean survival time, and 
 * quantiles). Either a vector of quantiles (i.e., times) or a vector of 
 * probabilities (when computing quantiles).  
 * @return An R data frame with the following columns:
 * - @c curve: An integer denoting a summarized curve. 
 * - @c sample: An integer denoting a randomly sampled parameter set.
 * - @c strategy_id: An integer denoting a treatemnt strategy id.
 * - @c patient_id: An integer denoting a patient id.
 * - @c x: The value of @p x.
 * - @c value: The summarized value (hazard, cumulative hazard, survival proportion,
 *     restricted mean survival time, or quantile).
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::DataFrame C_psm_curves_summary(Rcpp::Environment R_PsmCurves, 
                                     std::vector<double> x, 
                                     std::string type, double dr){
  // Initialize
  std::unique_ptr<hesim::psm::surv_mods> survmods = hesim::psm::surv_mods::create(R_PsmCurves);
  int n_models = survmods->get_n_models();
  int n_samples = survmods->get_n_samples();
  int n_obs = survmods->get_n_obs();
  int N = n_models * n_samples * n_obs * x.size();
  hesim::psm::surv_summary out(N);      
  
  // Loop
  int counter = 0;
  for (int m = 0; m < n_models; ++m){
    for (int s = 0; s < n_samples; ++s){
      for (int i = 0; i < n_obs; ++i){
        std::vector<double> res_vec;
        if (type == "survival"){
          res_vec = survmods->summary(m, s, i, x, "survival");
        }
        else if (type == "cumhazard"){
          res_vec = survmods->summary(m, s, i, x, "cumhazard");
        }
        else if (type == "hazard"){
          res_vec = survmods->summary(m, s, i, x, "hazard");
        }
        else if (type == "rmst"){
          res_vec = survmods->summary(m, s, i, x, "rmst", dr); 
        }
        else if (type == "quantile"){
          res_vec = survmods->quantile(m, s, i, x);
        } 
        else{
          Rcpp::stop("Selected type for summarizing survival distribution is not available.");
        }
        
        // store vector
        for (int j = 0; j < x.size(); ++j){
          out.curve_[counter] = m;
          out.sample_[counter] = s;
          out.strategy_id_[counter] = survmods->strategy_id_[i];
          out.patient_id_[counter] = survmods->patient_id_[i];
          out.x_[counter] = x[j];
          out.value_[counter] = res_vec[j];
          ++counter;
        } // end loop over health states
      } // end loop over observations
    } // end loop over samples
  } // end loop over models

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

/***************************************************************************//** 
 * @ingroup psm
 * Simulate health state probabilities with a partitioned survival model.
 * This function is exported to @c R and used in @c Psm$sim_stateprobs().
 * @param R_psm An R object of class @c Psm
 * @return The same output returned by hesim::psm::stateprobs::sim.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::List C_psm_sim_stateprobs(Rcpp::Environment R_Psm){
  hesim::psm::stateprobs stprobs(R_Psm);
  return stprobs.sim();
}

/***************************************************************************//** 
 * @ingroup psm
 * Simulate weighted length of stay from health state probabilities simulated
 * using a partitioned survival model.
 * This function is exported to @c R and used in @c Psm$sim_costs() and
 * @c Psm$sim_qalys(). 
 * @param R_psm An @c R object of class @c Psm
 * @param R_stateprobs State probabilities computed using an R object of class 
 * @c Psm. (This is needed in addition to @p R_psm because a modified copy of
 * @c Psm$stateprobs_ is passed to @c C++ and its more efficient to copy a single
 * data member than the entire class object.)
 * @param dr Discount rate.
 * @param type "costs" for costs; "qalys" for quality-adjusted life-years.
 * @param categories Categories with a given @p type. For QALYs, there is only one
 * category ("qalys"), but for costs this could consist of different cost categories
 * such as drug acquisition and administration costs, resource use costs, etc. 
 * @return An @c R data frame with columns equivalent to the data members in
 * hesim::wlos_out_out.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::DataFrame C_psm_sim_wlos(Rcpp::Environment R_Psm, Rcpp::DataFrame R_stateprobs, 
                              std::vector<double> dr, std::string type,
                              std::vector<std::string> categories){
  hesim::wlos wlos(R_Psm, type);
  hesim::stateprobs_out stprobs(R_stateprobs);
  std::vector<double> times = Rcpp::as<std::vector<double> > (R_Psm["t_"]);
  hesim::wlos_out out = wlos(stprobs, times, dr, categories);
  return out.create_R_data_frame();
}




