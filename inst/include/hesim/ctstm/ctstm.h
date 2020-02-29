# ifndef HESIM_CTSTM_H
# define HESIM_CTSTM_H
#include <hesim/statmods/obs_index.h>
#include <hesim/statmods/statmods.h>

namespace hesim {

/** @ingroup ctstm 
 * Classes and functions for simulating continuous time state transiton models.
 */
namespace ctstm {

/***************************************************************************//** 
 * Statistical models for health state transitions.
 ******************************************************************************/
class transmod{
public:
  statmods::obs_index obs_index_; ///< A statmods::obs_index object.
  trans_mat trans_mat_; ///<A transition matrix object. 
  
  /** 
   * The constructor.
   * Instantiates the base class. 
   * @param R_CtstmTrans An @c R object of class @c CtstmTrans.
   */ 
  transmod(Rcpp::Environment R_CtstmTrans)
    : obs_index_(Rcpp::as<Rcpp::List>(R_CtstmTrans["input_data"])),
      trans_mat_(Rcpp::as<arma::mat>(R_CtstmTrans["trans_mat"])){
  }
  virtual ~transmod() {}
  
  static std::unique_ptr<transmod> create(Rcpp::Environment R_CtstmTrans);
  
  /** 
   * Return the number of modeled treatment strategies.
   */    
  int get_n_strategies() const {
    return obs_index_.n_strategies_;
  }
  
  /** 
   * Return the number of modeled patients.
   */     
  int get_n_patients() const {
    return obs_index_.n_patients_;
  }
  
  /** 
   * Return the number of health states in the multi-state model.
   */    
  int get_n_states() const {
    return obs_index_.n_healthvals_;
  }
  
  /** 
   * Return the number of health state transitions.
   */     
  int get_n_transitions() {
    return trans_mat_.n_trans_;
  }
  
  /** 
   * Return the number of randomly sampled parmeter sets.
   */       
  int virtual get_n_samples() = 0;
  
  /** 
   * Set the maximum of the support for the survival model.
   */    
  void virtual set_max_x(double max_x) = 0;
  
  /** 
   * Summarize the survival model for a given health state transition.
   * The summary is conditional on the randomly sampled parameter set,
   * the health state transition, the treatment strategy, and patient. The
   * member obs_index_ should be updated by treatment strategy, and patient
   * before calling summary. 
   * @param trans A health state transition identification number.
   * @param sample A random sample of the parameters from the posterior
   * distribution.tes
   * @param t Times at which to make predictions.
   * @param type "hazard" for hazards; "cumhazard" for cumulative hazards.
   * @return Summary measure as determined by @p type for each time @p t.
   */    
  std::vector<double> virtual summary(int trans, int sample, std::vector<double> t, 
                                      std::string type) = 0;
  
  /** 
   * Randomly sample a health state transition.
   * Randomly sample an r->s health state transition; that is, a transition
   * from health state r to health state s.
   */    
  double virtual random(int trans, int sample) = 0;
  
  /** 
   * Randomly sample a health state transition from a truncated distributon.
   * Randomly sample an r->s health state transition; that is, a transition
   * from health state r to health state s that is left truncated at @c lower.
   */    
  double virtual trandom(int trans, int sample, double lower) = 0;  
};

/***************************************************************************//**
 * A multi-state model with a joint likelihood.
 ******************************************************************************/ 
class mstate : public transmod {
private:
  statmods::surv survmod_;
  
  /** 
   * Initialize survmod_ from an @c R object of class @c CtstmTrans. 
   */   
  static statmods::surv init_survmod_(Rcpp::Environment R_CtstmTrans){

    // Input matrices
    Rcpp::List R_input_mats = Rcpp::as<Rcpp::List> (R_CtstmTrans["input_data"]);
    Rcpp::List X_list = R_input_mats["X"];
    vecmats X = Rcpp::as<vecmats>(X_list);

    // Create statistical models
    Rcpp::List R_params_list = Rcpp::as<Rcpp::List>(R_CtstmTrans["params"]);
    statmods::params_surv  params_surv(R_params_list);
    statmods::surv survmod(X, params_surv);
    
    // Return
    return survmod;
  }
  
public:
 /** 
   * The constructor.
   * Instantiates the multi-state model.
   */ 
  mstate(Rcpp::Environment R_CtstmTrans) 
    : transmod(R_CtstmTrans),
    survmod_(init_survmod_(R_CtstmTrans)) {}

  int get_n_samples() {
    return survmod_.get_n_samples();  
  }
  
  void set_max_x(double max_x){
    survmod_.dist_->max_x_ = max_x;
  }
  
  std::vector<double> summary(int trans, int sample, std::vector<double> t, 
                              std::string type) {
    obs_index_.set_health_index(trans);
    return survmod_.summary(sample, obs_index_(), t, type);
  }
  
  double random(int trans, int sample) {
    obs_index_.set_health_index(trans);
    return survmod_.random(sample, obs_index_());
  }
  
  double trandom(int trans, int sample, double lower) {
    obs_index_.set_health_index(trans);
    return survmod_.trandom(sample, obs_index_(), lower, INFINITY);
  }  
    
};

/***************************************************************************//**
 * A multi-state model with transition specific models.
 ******************************************************************************/ 
class mstate_list : public transmod {
private:
  std::vector<statmods::surv> survmods_; ///< A vector of survival models for each transition.
  
  /** 
   * Initialize the survmods_ vector from an @c R object of class @c CtstmTrans. 
   */   
  static std::vector<statmods::surv> init_survmods_(Rcpp::Environment R_CtstmTrans){

    // Input matrices
    Rcpp::List R_input_mats = Rcpp::as<Rcpp::List> (R_CtstmTrans["input_data"]);
    Rcpp::List X_list = R_input_mats["X"];
    vecmats_2d X = Rcpp::as<vecmats_2d>(X_list);

    // Create statistical models
    Rcpp::List R_params_list = Rcpp::as<Rcpp::List>(R_CtstmTrans["params"]);
    std::vector<statmods::surv> survmods;
    for (int i = 0; i <R_params_list.size(); ++i){
      Rcpp::List R_params_list_i = Rcpp::as<Rcpp::List>(R_params_list[i]);
      statmods::params_surv  params_surv_i(R_params_list_i);
      statmods::surv survmod_i(X.at(i), params_surv_i);
      survmods.push_back(std::move(survmod_i));
    }
    return survmods;
  }
  
public:
 /** 
   * The constructor.
   * Instantiates the multi-state model.
   */ 
  mstate_list(Rcpp::Environment R_CtstmTrans) 
    : transmod(R_CtstmTrans),
    survmods_(init_survmods_(R_CtstmTrans)) {}
  
  int get_n_samples() {
    return survmods_[0].get_n_samples();  
  }
  
  void set_max_x(double max_x){
    int n_models = survmods_.size();
    for (int i = 0; i < n_models; ++i){
      survmods_[i].dist_->max_x_ = max_x;
    }
  }  
  
  std::vector<double> summary(int trans, int sample,
                              std::vector<double> t, std::string type) {
    return survmods_[trans].summary(sample, obs_index_(), t, type);
  }
  
  double random(int trans, int sample) {
    return survmods_[trans].random(sample, obs_index_());
  }
  
  double trandom(int trans, int sample, double lower) {
    return survmods_[trans].trandom(sample, obs_index_(), lower, INFINITY);
  } 
};

 /** 
   * A factory function.
   * Creates a multi-state model of health state transitions of the class specified
   * in @c R_CtstmTrans. The class must be inherited from ctstm::transmod.
   * @param R_CtstmTrans An @c R object of class @c CtstmTrans.
   * @return A unique pointer to the abstract base class ctstm::transmod. 
   */ 
inline std::unique_ptr<transmod> transmod::create(Rcpp::Environment R_CtstmTrans) {
  Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_CtstmTrans["params"]);
  transmod * mod;
  if (Rf_inherits(R_params, "params_surv_list")){
    mod = new mstate_list(R_CtstmTrans); 
  }
  else if (Rf_inherits(R_params, "params_surv")){
    mod = new mstate(R_CtstmTrans); 
  }
  else{
    Rcpp::stop("The selected statistical model is not available.");
  }
  std::unique_ptr<transmod> uptr(mod);
  return uptr;
}

/***************************************************************************//** 
 * Data container for storing summaries of models of health state transitions.
 ******************************************************************************/ 
struct transmod_summary{
  std::vector<int> trans_; ///< The health state transition.
  std::vector<int> sample_; ///< A randomly sampled parameter set.
  std::vector<int> strategy_id_; ///< The treatment strategy ID.
  std::vector<int> patient_id_; ///< The patient ID.
  std::vector<int> grp_id_; ///< The subgroup ID.
  std::vector<double> t_; ///< The time. 
  std::vector<double> value_; ///< The summarized value (hazard or cumulative hazard).
  
  /** 
   * A default constructor.
   * Instantiates a data container for a predicted survival curve.
   */ 
  transmod_summary() {};
  
/** 
   * A constructor.
   * Instantiates a data container for a predicted survival curve where all 
   * vectors in the container are initialized to a size @c n.
   */
  transmod_summary(int n) {
    trans_.resize(n);
    sample_.resize(n);
    strategy_id_.resize(n);
    patient_id_.resize(n);
    grp_id_.resize(n);
    t_.resize(n);
    value_.resize(n);
  }  
};

} // end namespace ctstm

} // end namespace hesim

# endif
