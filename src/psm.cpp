#include <hesim/psm.h>
#include <hesim/statevals.h>


namespace hesim {

namespace psm {

/****************
* Survival models
****************/
// Base case
surv_mods::surv_mods(Rcpp::Environment R_PsmCurves)
  : obs_index_(Rcpp::as<Rcpp::List>(R_PsmCurves["input_mats"])){
  Rcpp::Environment R_input_mats = Rcpp::as<Rcpp::Environment > (R_PsmCurves["input_mats"]);
  strategy_id_ = Rcpp::as<std::vector<int> >(R_input_mats["strategy_id"]);
  patient_id_ = Rcpp::as<std::vector<int> >(R_input_mats["patient_id"]);
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
    Rcpp::List R_input_mats = Rcpp::as<Rcpp::List > (R_PsmCurves["input_mats"]);
    X_ = Rcpp::as<vecmats_2d>(R_input_mats["X"]);
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
  int n_strategies = survmods->obs_index_.n_strategies_;
  int n_patients = survmods->obs_index_.n_patients_;
  int N = n_models * n_samples * n_strategies * n_patients * x.size();
  hesim::psm::surv_summary out(N);

  //Loop
  int counter = 0;
  for (int s = 0; s < n_samples; ++s) { // samples loop
    for (int k = 0; k < n_strategies; ++k) { // strategies loop
    survmods->obs_index_.set_strategy_index(k);
      for (int i = 0; i < n_patients; ++i) { // patients loop
        survmods->obs_index_.set_patient_index(i);
        for (int h = 0; h < n_models; ++h) { // survival models loop

          // Summarize each survival model
          int obs = survmods->obs_index_();
          std::vector<double> res_vec;
          if (type == "survival"){
            res_vec = survmods->summary(h, s, obs, x, "survival");
          }
          else if (type == "cumhazard"){
            res_vec = survmods->summary(h, s, obs, x, "cumhazard");
          }
          else if (type == "hazard"){
            res_vec = survmods->summary(h, s, obs, x, "hazard");
          }
          else if (type == "rmst"){
            res_vec = survmods->summary(h, s, obs, x, "rmst", dr);
          }
          else if (type == "quantile"){
            res_vec = survmods->quantile(h, s, obs, x);
          }
          else{
            Rcpp::stop("Selected type for summarizing survival distribution is not available.");
          }

          // Store results
          for (int j = 0; j < x.size(); ++j){
            out.curve_[counter] = h;
            out.sample_[counter] = s;
            out.strategy_id_[counter] = survmods->obs_index_.get_strategy_id();
            out.patient_id_[counter] = survmods->obs_index_.get_patient_id();
            out.x_[counter] = x[j];
            out.value_[counter] = res_vec[j];
            ++counter;
          }

        } // end survival models loop
      } // end patient loop
    } // end strategoes loop
  } // end samples loop

  // Return
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["sample"] = out.sample_,
    Rcpp::_["strategy_id"] = out.strategy_id_,
    Rcpp::_["patient_id"] = out.patient_id_,
    Rcpp::_["curve"] = out.curve_,
    Rcpp::_["x"] = out.x_,
    Rcpp::_["value"] = out.value_,
    Rcpp::_["stringsAsFactors"] = false
  );
  return out_df;
}

/***************************************************************************//** 
 * @ingroup psm
 * Helper function to compute helath state probabilities for a single row in
 * @c survival_ from the @c R class @c Psm. Used in C_psm_sim_stateprobs.
 * @param index An index denoting a row number in @c survival_.
 * @param state The health state to compute a probability for.
 * @param n_states The number of health states in the partitioned survival
 * model.
 * @param n_times The number of times at which survival was simulated. 
 * @param[out] n_crossings Incremented by 1 if the survival curves cross for this
 * row.
 * @param[in] survcurves_ The simulated survival curves passed from @c R to @c
 * C++.
 * @return The health state probability. 
 ******************************************************************************/ 
double stateprobs_sim1(int index, int state, int n_states, int n_times,
                       int &n_crossings, hesim::psm::surv_summary &survcurves) {
  if (state == 0){
    return survcurves.value_[index];
  }
  else if (state > 0 && state < (n_states - 1)){
    double survival = survcurves.value_[index];
    double survival_prior = survcurves.value_[index - n_times]; // Each health state is repeated n_times.
    if (survival < survival_prior){
      survival = survival_prior;
      ++n_crossings;
    }
    return survival - survival_prior;
  } 
  else{
    return 1 - survcurves.value_[index];
  }
}

/***************************************************************************//** 
 * @ingroup psm
 * Simulate health state probabilities with a partitioned survival model.
 * This function is exported to @c R and used in @c Psm$sim_stateprobs().
 * @param R_psm_survival The @c survival_ member of an R object of class @c Psm.
 * @param n_samples The number of random samples of the parameter set.
 * @param n_strategies The number of treatment strategies.
 * @param n_patients The number of patients.
 * @param n_states The number of health states.
 * @param n_times The number of times at which survival was simulated. 
 * @return The same output returned by hesim::psm::stateprobs::sim.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::List C_psm_sim_stateprobs(Rcpp::DataFrame R_psm_survival,
                                int n_samples, int n_strategies, int n_patients,
                                int n_states, int n_times){
  hesim::psm::surv_summary surv_curves(R_psm_survival); 
  hesim::stateprobs_out out(n_samples * n_strategies * n_patients * n_states * n_times);
  int n_curves = n_states - 1;
  
  // Loop
  int n_crossings = 0;
  int counter = 0;
  for (int s = 0; s < n_samples; ++s){ // begin samples loop
    for (int k = 0; k < n_strategies; ++k) { // begin strategies loop
      for (int i = 0; i < n_patients; ++i){ // begin patient loop
        for (int h = 0; h < n_states; ++h){ // begin loop over health states
          for (int t = 0; t < n_times; ++t) { // begin loop over time
            
            
          // Get the index for the survival curves
          int sc = 0;
          if (h < (n_states - 1)) { // If not in last health state
            sc = h;
          }
          else{
            sc = h - 1; // If the last health state
          }
          int index = s * n_strategies * n_patients * n_curves * n_times + 
                              k * n_patients * n_curves * n_times +
                              i * n_curves * n_times +
                              sc * n_times + 
                              t;
          
          // Results
          out.sample_[counter] = surv_curves.sample_[index];
          out.strategy_id_[counter] = surv_curves.strategy_id_[index];
          out.patient_id_[counter] = surv_curves.patient_id_[index];
          out.state_id_[counter] = h;
          out.t_[counter] = surv_curves.x_[index];
          out.prob_[counter] = stateprobs_sim1(index, h, n_states, n_times,
                                               n_crossings, surv_curves);          
            
          ++counter;
          }
        } // end health state loop
      } // end patient loop
    } // end treatment strategies loop
  } // end sample loop

  // Return
  Rcpp::DataFrame stateprobs_df = out.create_R_data_frame();
  
  return(Rcpp::List::create(
    Rcpp::_["stateprobs"] = stateprobs_df,
    Rcpp::_["n_crossings"] = n_crossings
  ));    
  
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
 * @c Psm$stateprobs_ is passed to @c C++ and it''s more efficient to copy a single
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




