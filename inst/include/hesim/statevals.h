# ifndef HESIM_STATEVALS_H
# define HESIM_STATEVALS_H
#include <hesim/statmods.h>
#include <hesim/input_data.h>
#include <hesim/integrate.h>

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
    Rcpp::List R_data = Rcpp::as<Rcpp::List>(R_StateVals["data"]);
    Rcpp::List X_list = Rcpp::as<Rcpp::List>(R_data["X"]);
    Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_StateVals["params"]);
    if (Rf_inherits(R_params, "params_lm")){
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
  /** 
   * The constructor.
   * Instantiates a model for simulating state values.
   * @param R_StateVals An @c R object from the @c hesim package of class
   * "StateVals".
   */ 
  statevals(Rcpp::Environment R_StateVals)
    : statmod_(init_statmod_(R_StateVals)) {}  

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
  std::vector<int> state_id_; ///< The health state id.
  std::vector<int> sample_; ///< The random sample of the parameters.
  std::vector<int> strategy_id_; ///< The treatment strategy id.
  std::vector<int> patient_id_; ////< The patient id.
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
    t_ = Rcpp::as<std::vector<double> >(R_stateprobs["t"]);
    prob_ = Rcpp::as<std::vector<double> >(R_stateprobs["prob"]);
  
    // R to C++ indexing
    add_constant(state_id_, -1);
    add_constant(sample_, -1);
    add_constant(strategy_id_, -1);
    add_constant(patient_id_, -1);
  }
  
  /** 
   * Create a data frame to pass to @c R.
   */   
  Rcpp::DataFrame create_R_data_frame(){
    return Rcpp::DataFrame::create(
      Rcpp::_["state_id"] = state_id_,
      Rcpp::_["sample"] = sample_,
      Rcpp::_["strategy_id"] = strategy_id_,
      Rcpp::_["patient_id"] = patient_id_,
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
struct wlos_out {
  std::vector<int> state_id_; ///< The health state id.
  std::vector<int> sample_; ///< The random sample of the parameters.
  std::vector<int> strategy_id_; ///< The treatment strategy id.
  std::vector<int> patient_id_; ////< The patient id.
  std::vector<double> dr_; ///< The discount rate.
  std::vector<std::string> category_; ///< For costs, the name of the cost category; 
                                     ///<for QALYs, simply 'qalys'.
  std::vector<double> value_; ///< The value of weighted length of stay.
  
  /** 
   * The constructor.
   * Instantiates a data container for storing output.
   * @param n Each member variable is initialized to have length @p n.
   */ 
  wlos_out(int n){
    state_id_.resize(n);
    sample_.resize(n);
    strategy_id_.resize(n);
    patient_id_.resize(n);
    dr_.resize(n);
    category_.resize(n);
    value_.resize(n);
  }
  
  /** 
   * Create a data frame to pass to @c R.
   */   
  Rcpp::DataFrame create_R_data_frame(){
    return Rcpp::DataFrame::create(
      Rcpp::_["state_id"] = state_id_,
      Rcpp::_["sample"] = sample_,
      Rcpp::_["strategy_id"] = strategy_id_,
      Rcpp::_["patient_id"] = patient_id_,
      Rcpp::_["dr"] = dr_,
      Rcpp::_["category"] = category_,
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
class wlos {
private:
  std::vector<statevals> statevals_; ///< A vector of models for simulating state values;
                                     ///< the vector will be of length 1 for 'qalys' and 
                                     ///< equal to the number of cost categories for 'costs'.
  statmods::obs_index obs_index_; ///< The @c obs_index class. 
  
  
  /** 
   * Initialize obs_index_.
   * Initialize obs_index_ from an @c R based @c hesim simulation model. 
   * @param R_model An @c R based @c hesim simulation model.
   */ 
  static statmods::obs_index init_obs_index_(Rcpp::Environment R_model, std::string type){
    if (type == "qalys"){
      Rcpp::Environment R_utility_model = Rcpp::as<Rcpp::Environment>(R_model["utility_model"]);
      return statmods::obs_index(Rcpp::as<Rcpp::List>(R_utility_model["data"]));
    }
    else if (type == "costs"){
      Rcpp::List R_cost_models = Rcpp::as<Rcpp::List>(R_model["cost_models"]);
      Rcpp::Environment R_cost_model_0 = Rcpp::as<Rcpp::Environment>(R_cost_models[0]);
      return statmods::obs_index(Rcpp::as<Rcpp::List>(R_cost_model_0["data"]));
    }
    else{
      Rcpp::stop("Values of 'costs' or 'qalys' can only be simulated.");
    }
  }

  /** 
   * Initialize statevals_.
   * Initialize statevals_ from an @c R based @c hesim simulation model. 
   * @param R_model An @c R based @c hesim simulation model.
   */   
  static std::vector<statevals> init_statevals_(Rcpp::Environment R_model, 
                                                std::string type){
    std::vector<statevals> statevals_vec;
    if (type == "qalys"){
      Rcpp::Environment R_utility_model = Rcpp::as<Rcpp::Environment>(R_model["utility_model"]);
      statevals stvals(R_utility_model);
      statevals_vec.push_back(std::move(stvals));
    } 
    else if (type == "costs"){
      Rcpp::List R_costs_models = Rcpp::as<Rcpp::List>((R_model["cost_models"]));
      for (int i = 0; i < R_costs_models.size(); ++i){
        Rcpp::Environment R_cost_model_i = Rcpp::as<Rcpp::Environment>(R_costs_models[i]);
        statevals stvals(R_cost_model_i);
        statevals_vec.push_back(std::move(stvals));
      }
    }
    else{
      Rcpp::stop("Values of 'costs' or 'qalys' can only be simulated.");
    }
    return statevals_vec;
  }
  
  /** 
   * Integrate probabilities.
   * Integrate health state probabilities weighted by the discount rate and 
   * assigned state values using the trapezoid rule (see @c math::trapz).
   * @param sample A random sample of the parameters from the posterior
   * distribution.
   * @param obs The observation (i.e., row index) for which to make a prediction
   * from the input matrix (or matrices when there are multiple parameters).
   * @param t_first, t_last The start and end times at which health state 
   * probabilities were simulated.
   * @param stateprob_first The beginning of health state probability values.
   * @param statevals A statistical model used to simulate the state values to assign
   * to health state. 
   * @return Weighted length of stay for a given parameter sample and observation.
   */ 
  double integrate_trapz(int sample, int obs,
                         std::vector<double>::iterator t_first, std::vector<double>::iterator t_last,
                         std::vector<double>::iterator stateprob_first, statevals &statevals,
                         double dr, std::string sim_type = "predict") {
    std::vector<double> value(std::distance(t_first, t_last));
    std::vector<double>::iterator stateprob_it = stateprob_first;
    for (std::vector<double>::iterator t_it = t_first; t_it != t_last; ++t_it){
      double statval = statevals.sim(sample, obs, sim_type);
      value[t_it - t_first] = exp(-dr * *t_it) * statval * *stateprob_it;
      ++stateprob_it;
    }
    return math::trapz(t_first, t_last, value.begin());
  }
                                                       
public:
  /** 
   * The constructor.
   * Instantiates an object for computing weighted length of stay.
   * @param R_model An @c R based @c hesim simulation model of class @c R6.
   * @param type 'costs' for costs or 'qalys' for QALYs.
   */ 
  wlos(Rcpp::Environment R_model, std::string type)
    : statevals_(init_statevals_(R_model, type)),
      obs_index_(init_obs_index_(R_model, type)){
  }
  
  /** 
   * Compute weighted length of stay.
   * Compute weighted length of stay given health state probabilities 
   * previously simulated at distinct times.
   * @param stateprobs
   * @param times
   * @param dr
   * @param categories
   * @return 
   */   
  wlos_out operator()(stateprobs_out stateprobs, 
                      std::vector<double> times,
                      std::vector<double> dr, 
                      std::vector<std::string> categories) {
    int N = stateprobs.prob_.size()/times.size();
    wlos_out out(N * dr.size() * statevals_.size());
    
    int counter = 0;
    for (int k = 0; k < statevals_.size(); ++k){
      for (int j = 0; j < dr.size(); ++j){
        int index = 0;
        double dr_j = dr[j];
        for (int i = 0; i < N; ++i){
          int state_id = stateprobs.state_id_[index];
          int strategy_id = stateprobs.strategy_id_[index];
          int patient_id = stateprobs.patient_id_[index];
          int sample = stateprobs.sample_[index];
          
          obs_index_.set_strategy_id(strategy_id);
          obs_index_.set_patient_id(patient_id);
          obs_index_.set_health_id(state_id);
          
          out.state_id_[counter] = state_id;
          out.sample_[counter] = sample;
          out.strategy_id_[counter] = strategy_id;
          out.patient_id_[counter] = patient_id;
          out.dr_[counter] = dr_j;
          out.category_[counter] = categories[k];
          out.value_[counter] = integrate_trapz(sample, obs_index_(),
                                                times.begin(), times.end(),
                                                stateprobs.prob_.begin() + index,
                                                statevals_[k],
                                                dr_j);
          index = index + times.size();
          ++counter;
        } // end loop over state probabilities
      } // end loop over discount rates
    } // end loop over state value models
    return out;
  }
  
}; // end class wlos


} // end namespace hesim


# endif
