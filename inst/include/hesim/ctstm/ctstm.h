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
 * Transition matrix.
 * A class for summarizing possible health state transitions in a multi-state 
 * model. 
 ******************************************************************************/
class trans_mat {
private:
  std::vector<std::vector<int> > trans_id_; ///< A vector of vectors. The outer vector 
                                            ///< denotes the starting state and each inner 
                                            ///< vector denotes the transition id
                                           ///< (indexed from 1 to patient::n_trans_)
                                           ///< corresponding to the possible transitions
                                           ///< from that state.
                                        
  std::vector<std::vector<int> > to_; ///< A vector of vectors. The outer vector
                                      ///< denotes the starting state and each inner
                                      ///< vector denotes a state that can be 
                                      ///< transitioned to.
                                        
  /** 
   * Count the number of non missing elements in the matrix. 
   * The number of non missing elements is equal to the number of 
   * possible transitions. 
   * @param m The same transition matrix as in the constructor.
   */                                          
  int count_non_nan(arma::mat m){
    int sum_non_nan = 0;
    for (int i = 0; i < m.n_rows; ++i){
      for (int j = 0; j < m.n_cols; ++j){
        if (!isnan(m(i, j))){
          ++sum_non_nan;
        }
      } // end loop over columns
    } // end loop over rows
    return sum_non_nan;
  }
  
  /** 
   * Determine whether each health state is absorbing.
   * @param trans A vector of vectors of the same format as trans_mat::trans_id_. Should
   * only be called after trans_ has been initialized. 
   */     
  std::vector<bool> is_absorbing(std::vector<std::vector<int> > trans){
    std::vector<bool> absorbing(trans.size());
    for (int i = 0; i < trans.size(); ++i){
      if(trans[i].size() > 0){
        absorbing[i] = false;
      }
      else {
        absorbing[i] = true;
      }
    } // end loop over states
    return absorbing;
  }
  
public:
  int n_trans_; ///< The total number of possible transitions.
  std::vector<bool> absorbing_; ///< A vector indicating whether each state is absorbing.
                               ///< A state is absorbing if a row in the transition matrix
                               ///< has all NAs. 
  
  /** 
   * The constructor.
   * @param m A matrix of integers indicating allowed transitions in a multi-state model
   *  in the format from the @c R package @c mstate. See
   * the argument "trans" in @c msprep in the @c mstate documentation.
   * @param R_index If TRUE, then transition ids in the matrix are assumed to be from R, 
   * and re-indexed to start from 0 (rather than 1).
   */                                       
  trans_mat(arma::mat m, bool R_index = true) {
    // Initialize n_trans_
    n_trans_ = count_non_nan(m);
    
    // Initialize trans_ and to_
    for (int i = 0; i < m.n_rows; ++i){
      arma::rowvec m_row = m.row(i);
      std::vector<int> trans_i;
      std::vector<int> to_i;
      for (int j = 0; j < m_row.n_elem; ++j){
        if(!isnan(m_row(j))){
          if (R_index){
            trans_i.push_back(m_row(j) - 1); 
          }
          else{
            trans_i.push_back(m_row(j)); 
          }
          to_i.push_back(j);
        }
      } // end loop over columns
      to_.push_back(to_i);
      trans_id_.push_back(trans_i);
    } // end loop over rows
    
    // Initialize absorbing_
    absorbing_ = is_absorbing(trans_id_);
  }
  
  /** 
   * Return transition number ids. 
   * @param from_state The state to transition from.
   * @reurn A vector of the transitions numbers from the specified health state.
   */   
  std::vector<int> trans_id(int from_state) {
    return trans_id_[from_state];
  }
  
  /** 
   * Return states that can be transitioned to. 
   * @param from_state The state to transition from.
   * @reurn A vector of the transition states that can be transitioned to 
   * from the specified health state.
   */   
  std::vector<int> to(int from_state) {
    return to_[from_state];
  }  
  
};

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
    : obs_index_(Rcpp::as<Rcpp::List>(R_CtstmTrans["data"])),
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
   * Return the number of modeled treatment lines.
   */    
  std::vector<int> get_n_lines() const {
    return obs_index_.n_lines_;
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
   * Summarize the survival model for a given health state transition.
   * The summary is conditional on the randomly sampled parameter set,
   * the health state transition, the treatment strategy, line, and patient. The
   * member obs_index_ should be updated by treatment strategy, line, and patient
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
    Rcpp::List data = Rcpp::as<Rcpp::List> (R_CtstmTrans["data"]);
    Rcpp::List X_list = data["X"];
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
  
  std::vector<double> summary(int trans, int sample, std::vector<double> t, 
                              std::string type) {
    obs_index_.set_health_index(trans);
    return survmod_.summary(sample, obs_index_(), t, type);
  }
  
  double random(int trans, int sample) {
    obs_index_.set_health_index(trans);
    return survmod_.random(sample, obs_index_());
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
    Rcpp::List data = Rcpp::as<Rcpp::List> (R_CtstmTrans["data"]);
    Rcpp::List X_list = data["X"];
    vecmats_2d X = Rcpp::as<vecmats_2d>(X_list);

    // Create statistical models
    Rcpp::List R_params_list = Rcpp::as<Rcpp::List>(R_CtstmTrans["params"]);
    std::vector<statmods::surv> survmods;
    for (int i = 0; i <R_params_list.size(); ++i){
      Rcpp::List R_params_list_i = Rcpp::as<Rcpp::List>(R_params_list[i]);
      statmods::params_surv  params_surv_i(R_params_list_i);
      statmods::surv survmod_i(X.at(0), params_surv_i);
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
  
  std::vector<double> summary(int trans, int sample,
                              std::vector<double> t, std::string type) {
    return survmods_[trans].summary(sample, obs_index_(), t, type);
  }
  
  double random(int trans, int sample) {
    return survmods_[trans].random(sample, obs_index_());
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
  std::vector<int> strategy_id_; ///< A treatment strategy id.
  std::vector<int> patient_id_; ///< A patient id.
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
    t_.resize(n);
    value_.resize(n);
  }  
};

/***************************************************************************//** 
 * Data container for storing health state probabilities.
 ******************************************************************************/ 
struct stateprobs_out{
  std::vector<int> sample_; ///< A randomly sampled parameter set.
  std::vector<int> strategy_id_; ///< A treatment strategy id.
  std::vector<double> state_id_; ///< The health state id. 
  std::vector<double> patient_id_; ///< The health state id. 
  std::vector<double> t_; ///< The time. 
  std::vector<double> prob_; ///< The health state probability. 
  
  /** 
   * A default constructor.
   * Instantiates a data container for storing simulated health state probabilities.
   */ 
  stateprobs_out() {};
  
/** 
   * A constructor.
   * Instantiates a data container for a predicted survival curve where all 
   * vectors in the container are initialized to a size @c n.
   */
  stateprobs_out(int n) {
    sample_.resize(n);
    strategy_id_.resize(n);
    state_id_.resize(n);
    patient_id_.resize(n);
    t_.resize(n);
    prob_.resize(n);
  }    
  
};


} // end namespace ctstm

} // end namespace hesim

# endif
