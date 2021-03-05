#include <hesim/psm.h>
#include <hesim/statevals.h>


namespace hesim {

namespace psm {

/****************
* Survival models
****************/
// Base case
surv_mods::surv_mods(Rcpp::Environment R_PsmCurves)
  : obs_index_(Rcpp::as<Rcpp::List>(R_PsmCurves["input_data"])){
  Rcpp::Environment R_input_mats = Rcpp::as<Rcpp::Environment > (R_PsmCurves["input_data"]);
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
    Rcpp::List R_input_mats = Rcpp::as<Rcpp::List > (R_PsmCurves["input_data"]);
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

class stateprobs{
private: 
  arma::cube surv_; ///< Survival curves
  int n_curves_; ///< The number of survival curves
  int n_states_; ///< Number of health states.
  int n_obs_; ///< Number of observations.
  int n_times_; ////< Number of times.

public:
  arma::cube prob_; ///< State probabilities
  std::vector<int> cross_; ///< Integers denoting whether survival curves crosses  
  
  /** 
   * The constructor.
   */ 
  stateprobs(arma::cube surv){
    surv_ = surv;
    n_obs_ = surv.n_slices;
    n_times_ = surv.n_rows;
    n_curves_ = surv.n_cols;
    n_states_ = n_curves_ + 1;
    cross_.resize(n_obs_, 0);
    prob_.resize(n_times_, n_states_, n_obs_);
  }
  
  /** 
   * Compute simulated state probabilities from survival curves for a single
   * observation and a single time point.
   * @param i An integer denoting the "observation the survival curve corresponds
   * to.
   * @param t An integer denoting the time at which a survival probability
   * was simulated.
   * @return None; @p prob_ and @p cross_ are updated.
   */   
  void sim(int i, int t) {
    
    // Probability in health state 1 is S1(t)
    prob_(t, 0, i) = surv_(t, 0, i);
    
    // Probability in health state 2, ... N-1 is Sn(t) − Sn−1(t)
    double surv_prior = surv_(t, 0, i);
    for (int j = 1; j < n_curves_; ++j) {
      if (surv_(t, j, i) >= surv_prior) { // Case where survival curves don't cross
        prob_(t, j, i) = surv_(t, j, i) - surv_prior;
        surv_prior = surv_(t, j, i);
      } else { // Case where survival curves cross
        prob_(t, j, i) = 0;
        ++cross_[i];
      }
    }

    // Probability in health state N is 1 - SN-1(t)
    prob_(t, n_curves_, i) = 1 - surv_prior;
  }
  
  /** 
   * Compute simulated state probabilities from survival curves across all 
   * observations and time points.
   * @return None; @p prob_ and @p cross_ are updated.
   */  
  void sim() {
    for (int i = 0; i < n_obs_; ++i) { // Loop over observations
      for (int t = 0; t < n_times_; ++t) { // Loop over times
        sim(i, t);
      }
    } // End loop over observations
  }
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
            out.grp_id_[counter] = survmods->obs_index_.get_grp_id();
            out.patient_wt_[counter] = survmods->obs_index_.get_patient_wt();
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
    Rcpp::_["grp_id"] = out.grp_id_,
    Rcpp::_["patient_wt"] = out.patient_wt_,
    Rcpp::_["curve"] = out.curve_,
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
 * @param surv A cube containing survival curves where rows are times, columns
 * denote the curve number, and slices are an observation (sorted by parameter sample,
 * treatment strategy, and patient).
 * @return A list with two elements. The first element @c prob contains health
 * state probabilities computed using partitioned survival analysis. It is of the
 * same dimensions as @c surv but has one more column (where columns denote the
 * health state rather than the curve number). The second element is an integer
 * valued vector with length equal to the number of slices in @c surv denoting
 * the number of times survival curves crossed for a given observation.
 ******************************************************************************/
// [[Rcpp::export]]
Rcpp::List C_psm_sim_stateprobs(arma::cube surv) {
  hesim::psm::stateprobs stprobs(surv);
  stprobs.sim();
  return(Rcpp::List::create(
    Rcpp::_["prob"] = stprobs.prob_,
    Rcpp::_["cross"] = stprobs.cross_
  )); 
}