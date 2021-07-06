# ifndef HESIM_STATEVALS_H
# define HESIM_STATEVALS_H

#include <hesim/statmods/statmods.h>
#include <hesim/statmods/obs_index.h>
#include <hesim/math/composite.h>

namespace hesim {

/***************************************************************************//** 
 * @ingroup general
 * Simulate health state values.
 * Simulate the values of health state by either predicting mean values
 * or sampling random values with a statistical model.
 ******************************************************************************/ 
class statevals {
private:
  /** 
   * Initialize @c statmod_.
   * Initialize the statistical model.
   * @param R_StateVals An object from the @c R package @c hesim of class "StateVals".
   * @return A unique pointer to statmods::statmod.
   */   
  static std::unique_ptr<statmods::statmod> init_statmod_(Rcpp::Environment R_StateVals){
    statmods::statmod *mod;
    Rcpp::List R_data;
    Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_StateVals["params"]);
    if (Rf_inherits(R_params, "tparams_mean")){
      mod = new statmods::pred_means(R_params);
    }
    else if (Rf_inherits(R_params, "params_lm")){
      R_data = Rcpp::as<Rcpp::List>(R_StateVals["input_data"]);
      Rcpp::List X_list = Rcpp::as<Rcpp::List>(R_data["X"]);      
      arma::mat X = Rcpp::as<arma::mat>(X_list["mu"]);
      statmods::params_lm params(R_params);
      mod = new statmods::lm(X, params);
    }
    else{
      Rcpp::stop("The class of 'params' is not supported.");
    }
    std::unique_ptr<statmods::statmod> uptr(mod);
    return uptr;
  }
  

                                            
public:
  std::unique_ptr<statmods::statmod> statmod_; ///< The statistical model used
                                            ///< to simulate state values.
  std::string method_; ///<The method used to sim ulate costs and quality-adjusted 
                       ///< life-years. See the @c R StateVals class.
                                            
  /** 
   * The constructor.
   * Instantiates a model for simulating state values.
   * @param R_StateVals An @c R object from the @c hesim package of class
   * "StateVals".
   */ 
  statevals(Rcpp::Environment R_StateVals)
    : statmod_(init_statmod_(R_StateVals)) {
    method_ = Rcpp::as<std::string>(R_StateVals["method"]);
  }  

  /** 
   * Simulate value of a health state.
   * @param sample A random sample of the parameters from the posterior
   * distribution.
   * @param obs The observation (i.e., row index) for which to make a prediction
   * from the input matrix (or matrices when there are multiple parameters).
   * @param type "predict" to predict mean values; "random" to randomly draw
   * a value from the probability distribution of the statistical model.
   * @return The simulated value.
   */ 
   double sim(int sample, int obs, std::string type = "predict") {
     if (type == "predict"){
       return statmod_->predict(sample,  obs);
     }
     else if (type == "random"){
       return statmod_->random(sample, obs);
     }
     else{
       Rcpp::stop("'type' must either be 'predict' or 'random'.");
     }
   }
  
}; // end statevals class

/***************************************************************************//** 
 * @ingroup general
 * Health state probability output.
 * Contains output of health state probability computations.
 ******************************************************************************/ 
struct stateprobs_out{
  std::vector<int> state_id_; ///< The health state ID.
  std::vector<int> sample_; ///< The random sample of the parameters.
  std::vector<int> strategy_id_; ///< The treatment strategy ID.
  std::vector<int> patient_id_; ////< The patient ID.
  std::vector<int> grp_id_; ///< The subgroup ID.
  std::vector<double> patient_wt_; ///< Weights given to patients.
  std::vector<double> t_; ///< Time.
  std::vector<double> prob_; ///< A health state probability.
  
  /** 
   * A default constructor.
   * Instantiates a data container where each member variable is an empty 
   * container with no elements.
   */ 
  stateprobs_out() {
  }
  
  /** 
   * A constructor.
   * Instantiates an empty data container for storing output.
   * @param n Each member variable is initialized to have length @p n.
   */ 
  stateprobs_out(int n) {
    state_id_.resize(n);
    sample_.resize(n);
    strategy_id_.resize(n);
    patient_id_.resize(n);
    grp_id_.resize(n);
    patient_wt_.resize(n);
    t_.resize(n);
    prob_.resize(n);    
  }
  
  /** 
   * A constructor.
   * Instantiates a filled data container with output passed from @c R.
   * @param R_stateprobs An @c R data frame with a column for each member variable
   * in the struct.
   */ 
  stateprobs_out(Rcpp::DataFrame R_stateprobs) {
    state_id_ = Rcpp::as<std::vector<int> >(R_stateprobs["state_id"]);
    sample_ = Rcpp::as<std::vector<int> >(R_stateprobs["sample"]);
    strategy_id_ = Rcpp::as<std::vector<int> >(R_stateprobs["strategy_id"]);
    patient_id_ = Rcpp::as<std::vector<int> >(R_stateprobs["patient_id"]);
    grp_id_ = Rcpp::as<std::vector<int> >(R_stateprobs["grp_id"]);
    if (!hesim::is_null(R_stateprobs, "patient_wt")){
      patient_wt_ = Rcpp::as<std::vector<double> >(R_stateprobs["patient_wt"]);
    } 
    else{
      patient_wt_.resize(patient_id_.size(), 1.0);
    }
    t_ = Rcpp::as<std::vector<double> >(R_stateprobs["t"]);
    prob_ = Rcpp::as<std::vector<double> >(R_stateprobs["prob"]);
  
    // R to C++ indexing
    add_constant(state_id_, -1);
    add_constant(sample_, -1);
    add_constant(strategy_id_, -1);
    add_constant(patient_id_, -1);
    add_constant(grp_id_, -1);
  }
  
  /** 
   * Create a data frame to pass to @c R.
   */   
  Rcpp::DataFrame create_R_data_frame(){
    return Rcpp::DataFrame::create(
      Rcpp::_["sample"] = sample_,
      Rcpp::_["strategy_id"] = strategy_id_,
      Rcpp::_["patient_id"] = patient_id_,
      Rcpp::_["grp_id"] = grp_id_,
      Rcpp::_["patient_wt"] = patient_wt_,
      Rcpp::_["state_id"] = state_id_,
      Rcpp::_["t"] = t_,
      Rcpp::_["prob"] = prob_,
      Rcpp::_["stringsAsFactors"] = false
    );
  }
  
}; // end struct stateprobs_out

/***************************************************************************//** 
 * @ingroup general
 * Weighted health state length of stay output.
 * Contains output of weighted health state length of stay computations.
 ******************************************************************************/ 
struct ev_out {
  std::vector<int> state_id_; ///< The health state ID.
  std::vector<int> sample_; ///< The random sample of the parameters.
  std::vector<int> strategy_id_; ///< The treatment strategy ID.
  std::vector<int> patient_id_; ////< The patient ID.
  std::vector<int> grp_id_; ////< The subgroup ID.
  std::vector<double> patient_wt_; ///< Weights given to patients.
  std::vector<double> dr_; ///< The discount rate.
  std::vector<std::string> outcome_; ///< The name of the outcome. 
                                     ///<For costs, the name of the cost category; 
                                     ///<for QALYs, simply 'qalys'.
  std::vector<double> value_; ///< The value of weighted length of stay.
  
  /** 
   * The constructor.
   * Instantiates a data container for storing output.
   * @param n Each member variable is initialized to have length @p n.
   */ 
  ev_out(int n){
    state_id_.resize(n);
    sample_.resize(n);
    strategy_id_.resize(n);
    patient_id_.resize(n);
    grp_id_.resize(n);
    patient_wt_.resize(n);
    dr_.resize(n);
    outcome_.resize(n);
    value_.resize(n);
  }
  
  /** 
   * Create a data frame to pass to @c R.
   */   
  Rcpp::DataFrame create_R_data_frame(){
    return Rcpp::DataFrame::create(
      Rcpp::_["sample"] = sample_,
      Rcpp::_["strategy_id"] = strategy_id_,
      Rcpp::_["patient_id"] = patient_id_,
      Rcpp::_["grp_id"] = grp_id_,
      Rcpp::_["patient_wt"] = patient_wt_,
      Rcpp::_["state_id"] = state_id_,
      Rcpp::_["dr"] = dr_,
      Rcpp::_["outcome"] = outcome_,
      Rcpp::_["value"] = value_,
      Rcpp::_["stringsAsFactors"] = false
    );
  }
  
}; // end class weighted_los_out 
  
/***************************************************************************//** 
 * @ingroup general
 * Compute weighted length of stay in health states.
 * Compute the weighted length of stay in each health state from time @c 0 to 
 * time @c T, which is defined as the integral of the probability of being in
 * that health state weighted by a discount rate @c dr and a state value assigned
 * to each health state from time @c0 to time @c T.
 ******************************************************************************/ 
class ev {
private:
  std::vector<statevals> statevals_; ///< A vector of models for simulating state values;
  
  /** 
   * Initialize obs_index_.
   * Initialize obs_index_ from an @c R based @c hesim simulation model. 
   * @param R_model An @c R based @c hesim simulation model.
   */ 
  static statmods::obs_index init_obs_index_(Rcpp::List R_statevals){
    Rcpp::Environment R_statevals_0 = Rcpp::as<Rcpp::Environment>(R_statevals[0]);
    return statmods::obs_index(Rcpp::as<Rcpp::List>(hesim::statmods::get_id_object(R_statevals[0])));    
  }  
  
  /** 
   * Initialize obs_indices_.
   * Initialize vector of @c obs_index_ objects for each outcome predicted
   * with an @c R based @c hesim simulation model. 
   * @param R_statevals A list of @c R based @c  models for state values by
   * outcome.
   */ 
  static std::vector<statmods::obs_index> init_obs_indices_(Rcpp::List R_statevals){
    std::vector<statmods::obs_index> obs_index_vec;
    for (int i = 0; i < R_statevals.size(); ++i){ 
      Rcpp::Environment R_statevals_i = Rcpp::as<Rcpp::Environment>(R_statevals[i]);
      obs_index_vec.push_back(statmods::obs_index(hesim::statmods::get_id_object(R_statevals_i)));
    }
    return obs_index_vec;
  }    
  
  /** 
   * Initialize statevals_.
   * Initialize statevals_ from an @c R based @c hesim state value model.  
   * @param R_statevals An @c R based @c hesim state value model.
   */   
  static std::vector<statevals> init_statevals_(Rcpp::List R_statevals){
    std::vector<statevals> statevals_vec;
    for (int i = 0; i < R_statevals.size(); ++i){
      Rcpp::Environment R_statevals_i = Rcpp::as<Rcpp::Environment>(R_statevals[i]);
      statevals stvals(R_statevals_i);
      statevals_vec.push_back(std::move(stvals));
      }
    return statevals_vec;
  }  
  
  /** 
   * Integrate given method.
   * Integrate over values at fixed time periods given chosen method.
   * @param times The fixed time period.
   * @param values_first A pointer to the beginning of the vector of values
   * associated with each time period.
   * @param method The method used to integrate state values. 
   */ 
   static double integrate(std::vector<double> &times, 
                    std::vector<double>::iterator values_first,
                    std::string method) {
     if (method == "trapz"){
       return math::trapz(times.begin(), times.end(), values_first);
     }
     else if (method == "riemann_left"){
       return math::riemann_left(times.begin(), times.end(), values_first);
     }
     else if (method == "riemann_right"){
       return math::riemann_right(times.begin(), times.end(), values_first);
     }
     else{
       Rcpp::stop("The selected integration method is not available.");
     }
   }
  
  /** 
   * Weighted length of stay.
   * Integrate health state probabilities weighted by the discount rate and 
   * assigned state values using the chosen method. 
   * Predictions are made for observations (i.e., row indices) determined
   * by @p obs_index. 
   * @param sample A random sample of the parameters from the posterior
   * distribution.
   * @param obs_index An object of class @c obs_index denoting the observation
   * index.
   * @param times Times at which state probabilities were computed. 
   * @param stateprob_first The beginning of health state probability values.
   * @param statevals A statistical model used to simulate the state values to assign
   * to health state. 
   * @param dr The discount rate. 
   * @param method The method used to integrate state values. 
   * @param sim_type "predict" for mean values or "random" for random samples.
   * @return Weighted length of stay for a given parameter sample and observation.
   */ 
  double sim_wlos(int sample, statmods::obs_index obs_index,
                  std::vector<double> times,
                  std::vector<double>::iterator stateprob_first, 
                  statevals &statevals,
                  double dr, std::string method = "trapz",
                  std::string sim_type = "predict") {
    
    // State values
    std::vector<double> value(std::distance(times.begin(), times.end()));
    auto stateprob_it = stateprob_first;
    auto t_start = times.begin();
    for (int t = 0 ; t < obs_index.n_times_; ++t){
      obs_index.set_time_index(t);
      auto t_stop = hesim::max_lt(t_start, times.end(), obs_index.get_time_stop());
      double statval = statevals.sim(sample, obs_index(), sim_type);
      for (auto t_it = t_start; t_it <= t_stop; ++t_it){
        value[t_it - times.begin()] = exp(-dr * *t_it) * statval * *stateprob_it;
        ++stateprob_it;
      } // Loop over times within time interval
      t_start = t_stop + 1;
    } // Loop over time intervals
    
    // Integrate
    return integrate(times, value.begin(), method);
  }

  /** 
   * One-time expected values.
   * Simulate expected values in a cohort model based on one-time state values
   * that occur when a patient enters a health state. Values only accrue at 
   * time 0.
   * @param sample A random sample of the parameters from the posterior
   * distribution.
   * @param obs_index An object of class @c obs_index denoting the observation
   * index.
   * @param stateprob_first The beginning of health state probability values.
   * @param statevals A statistical model used to simulate the state values to assign
   * to health state. 
   * @param sim_type "predict" for mean values or "random" for random samples.
   * @return Weighted length of stay for a given parameter sample and observation.
   */   
  double sim_starting(int sample, statmods::obs_index obs_index,
                      double &stateprob_start,
                      statevals &statevals,
                      std::string sim_type = "predict") {
    obs_index.set_time_index(0);
    double stval = statevals.sim(sample, obs_index(), sim_type);
    return stval * stateprob_start;
  }
                                                       
public:
  statmods::obs_index obs_index_; ///< The @c obs_index class. 
  std::vector<statmods::obs_index> obs_indices_; ///< The @c obs_index class.
  
  /** 
   * The constructor.
   * Instantiates an object for computing weighted length of stay.
   * @param R_model An @c R based @c hesim simulation model of class @c R6.
   * @param type 'costs' for costs or 'qalys' for QALYs.
   */ 
  ev(Rcpp::List R_statevals)
    : statevals_(init_statevals_(R_statevals)),
      obs_index_(init_obs_index_(R_statevals)),
      obs_indices_(init_obs_indices_(R_statevals)){
  }
  
  /** 
   * Length of stay.
   * Integrate health state probabilities (potentially weighted by the discount 
   * rate) using the chosen method.
   * @param times Times at which state probabilities were computed. 
   * @param stateprob_first The beginning of health state probability values.
   * @param dr The discount rate. 
   * @param method The method used to integrate state values. 
   * @return Length of stay associated with a vector of health state 
   * probabilities.
   */ 
  static double sim_los(std::vector<double> &times,
                        std::vector<double>::iterator stateprob_first, 
                        double dr, 
                        std::string method = "trapz") {
    
    auto stateprob_it = stateprob_first;
    std::vector<double> value(std::distance(times.begin(), times.end()));
    for (int i = 0; i < times.size(); ++i) {
      value[i] = exp(-dr * times[i]) * *stateprob_it;
      ++stateprob_it;
    }
    
    return integrate(times, value.begin(), method);
  }
  
  /** 
   * Simulate expected values.
   * Simulate costs and quality-adjusted life-years (QALYs) as a function
   * of simulated state occupancy probabilities.
   * previously simulated at distinct times.
   * @param stateprobs
   * @param times
   * @param dr
   * @param outcomes
   * @return 
   */   
  ev_out operator()(stateprobs_out stateprobs, 
                      std::vector<double> times,
                      std::vector<double> dr, 
                      std::vector<std::string> outcomes,
                      std::string integrate_method = "trapz") {
    int N = stateprobs.prob_.size()/times.size();
    int n_samples = statevals_[0].statmod_->get_n_samples();
    ev_out out(N * dr.size() * statevals_.size());
    
    int counter = 0;
    for (int k = 0; k < statevals_.size(); ++k){ // start outcome loop
      for (int j = 0; j < dr.size(); ++j){ // start discount rate loop
        int integrate_start = 0;
        double dr_j = dr[j];
        for (int s = 0; s < n_samples; ++s){ // start samples loop
          for (int ts = 0; ts < obs_indices_[k].n_strategies_; ++ts){ // start treatment strategies loop
          obs_indices_[k].set_strategy_index(ts);
            for (int p = 0; p < obs_indices_[k].n_patients_; ++ p){ // start patient loop
              obs_indices_[k].set_patient_index(p);
              for (int h = 0; h < obs_indices_[k].n_healthvals_; ++h){ // health state state loop
                obs_indices_[k].set_health_index(h);
                
                out.sample_[counter] = s;
                out.strategy_id_[counter] = obs_indices_[k].get_strategy_id();
                out.patient_id_[counter] = obs_indices_[k].get_patient_id();
                out.grp_id_[counter] = obs_indices_[k].get_grp_id();
                out.patient_wt_[counter] = obs_indices_[k].get_patient_wt();
                out.state_id_[counter] = obs_indices_[k].get_health_id();;
                out.dr_[counter] = dr_j;
                out.outcome_[counter] = outcomes[k];
                if (statevals_[k].method_ == "wlos"){
                  out.value_[counter] = sim_wlos(s, obs_indices_[k],
                                                    times,
                                                    stateprobs.prob_.begin() + integrate_start,
                                                    statevals_[k],
                                                    dr_j, integrate_method);
                } 
                else{
                  out.value_[counter] = sim_starting(s, obs_indices_[k],
                                                     stateprobs.prob_[integrate_start],
                                                     statevals_[k]);
                }
                integrate_start = integrate_start + times.size();
                ++counter;
              } // end health state loop
            } // end patient loop
          } // end strategy loop
        } // end samples loop
      } // end loop over discount rates
    } // end loop over outcomes
    return out;
  }
  
}; // end class ev


} // end namespace hesim


# endif
