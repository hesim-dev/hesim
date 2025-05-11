# ifndef HESIM_DTSTM_H
# define HESIM_DTSTM_H
#include <hesim/statmods/obs_index.h>
#include <hesim/statmods/statmods.h>

namespace hesim {

/** @ingroup dtstm 
 * Classes and functions for simulating discrete time state transiton models.
 */
namespace dtstm {

/***************************************************************************//** 
 * Transition probability base class
 ******************************************************************************/

class trans_model{
public:
  statmods::obs_index obs_index_; ///< A statmods::obs_index object.
  arma::rowvec start_stateprobs_;///< A vector with length equal to the number of 
  ///< health states containing the probability that 
  ////< the cohort is in each health state at the start of the simulation.
  
  /** 
   * The constructor.
   * Instantiates the base class. 
   * @param R_DtstmTrans An @c R object of class @c DtstmTrans.
   */ 
  trans_model(Rcpp::Environment R_CohortDtstmTrans)
    : obs_index_(hesim::statmods::get_id_object(R_CohortDtstmTrans)),
      start_stateprobs_(Rcpp::as<arma::rowvec>(R_CohortDtstmTrans["start_stateprobs"])){}
  virtual ~trans_model() {}  
  
  static std::unique_ptr<trans_model> create(Rcpp::Environment R_CohortDtstmTrans);
  
  virtual int get_n_states() = 0;
  
  virtual int get_n_samples() = 0;
  
  /** 
   * Transition probability matrix.
   * Return a transition probability matrix as a function of a parameter
   * sample and time.
   * @param s The parameter sample.
   * @param time The time.
   */     
  virtual arma::mat tpmatrix(int s, double time) = 0;
  
  void set_time_index(double time){
    while (time > obs_index_.get_time_stop()){
      obs_index_.set_time_index(obs_index_.get_time_index() + 1);
    }
  };
  
}; // end class trans_model

/***************************************************************************//**
 * Transition probabilities from a tparams_transprobs @c R object
 ******************************************************************************/ 
class tparams_transprobs : public trans_model {
private:
  int n_samples_; ///< The number of randomly sampled parameter sets.
  int n_states_; ///< The number of health states.
  int n_obs_; ///< The number of treatment strategy, patient, and time interval
  ///< observations.
  arma::cube value_; ///< An array of transition probability matrices.
  
public:
  /** 
   * The constructor.
   * Instantiates the object.
   */ 
  tparams_transprobs(Rcpp::Environment R_CohortDtstmTrans) 
    : trans_model(R_CohortDtstmTrans) {
    Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_CohortDtstmTrans["params"]);
    value_ = Rcpp::as<arma::cube>(R_params["value"]);
    n_samples_ = Rcpp::as<int>(R_params["n_samples"]);
    n_states_ = value_.slice(0).n_rows;
    n_obs_ = value_.n_slices/n_samples_;
  }
  
  int get_n_states(){
    return n_states_;
  }
  
  int get_n_samples(){
    return n_samples_;
  }
  
  arma::mat tpmatrix(int s, double time){
    set_time_index(time);
    return value_.slice(obs_index_() + s * n_obs_);
  }
  
}; // end class tparams_transprobs

/***************************************************************************//**
 * Transition probabilities from a params_mlogit_list @c R object
 ******************************************************************************/ 
class mlogit_list : public trans_model {
private:
  int n_states_; ///< The number of health states.
  std::vector<statmods::mlogit> mlogit_list_; ///< A vector of multinomial logit models.
  std::vector<bool> absorbing_; ///< Boolean for whether a health state is an absorbing state.
  arma::mat p_; ///< The transition matrix.
  int current_obs_index_;
  
  /** 
   * Initialize the mlogit_list_ vector from an @c R object of class @c CohortDtstmTrans. 
   */   
  static std::vector<statmods::mlogit> init_mlogit_list_(Rcpp::Environment R_CohortDtstmTrans){
    
    // Input matrices
    Rcpp::List R_input_mats = Rcpp::as<Rcpp::List> (R_CohortDtstmTrans["input_data"]);
    Rcpp::List X_list = R_input_mats["X"];
    vecmats X = Rcpp::as<vecmats>(X_list);
    
    // Transition matrix
    arma::mat trans_mat = Rcpp::as<arma::mat> (R_CohortDtstmTrans["trans_mat"]) - 1;
    
    // Create statistical models
    Rcpp::List R_params_list = Rcpp::as<Rcpp::List>(R_CohortDtstmTrans["params"]);
    std::vector<statmods::mlogit> mlogit_list;
    for (int i = 0; i < R_params_list.size(); ++i){
      Rcpp::List R_params_list_i = Rcpp::as<Rcpp::List>(R_params_list[i]);
      statmods::params_mlogit  params_mlogit_i(R_params_list_i);
      statmods::mlogit mlogit_list_i(X.at(i), params_mlogit_i, trans_mat.row(i));
      mlogit_list.push_back(std::move(mlogit_list_i));
    }
    return mlogit_list;
  }  
  
public:
  /** 
   * The constructor.
   * Instantiates the object.
   */ 
  mlogit_list(Rcpp::Environment R_CohortDtstmTrans) 
    : trans_model(R_CohortDtstmTrans),
      mlogit_list_(init_mlogit_list_(R_CohortDtstmTrans)) {
    trans_mat tmat(Rcpp::as<arma::mat> (R_CohortDtstmTrans["trans_mat"])); 
    absorbing_ = tmat.absorbing_;
    n_states_ = tmat.n_states_;
    p_.eye(n_states_, n_states_);
    current_obs_index_ = obs_index_();
  }
  
  int get_n_states(){
    return n_states_;
  }
  
  int get_n_samples(){
    return mlogit_list_.at(0).get_n_samples();
  }
  
  arma::mat tpmatrix(int s, double time){
    set_time_index(time);
    if (obs_index_() != current_obs_index_ || obs_index_() == 0){
      current_obs_index_ = obs_index_();
      for (int i = 0; i < (int) mlogit_list_.size(); ++i){
        if (!absorbing_[i]){
          p_.row(i) =  mlogit_list_[i].multi_predict(s, obs_index_());
        }
      }
    }
    return p_;
  }
  
}; // end class tparams_transprobs

/***************************************************************************//**
 * A factory function.
 * Creates a discrete time cohort transition model of the class specified
 * in @c R_CohortDtstmTrans. The class must be inherited from dtstm::trans_model.
 * @param R_CohortDtstmTrans An @c R object of class @c CohortDtstmTrans.
 * @return A unique pointer to the abstract base class dtstm::trans_model. 
 ******************************************************************************/ 
inline std::unique_ptr<trans_model> trans_model::create(Rcpp::Environment R_CohortDtstmTrans) {
  Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_CohortDtstmTrans["params"]);
  trans_model * mod;
  if (Rf_inherits(R_params, "tparams_transprobs")){
    mod = new tparams_transprobs(R_CohortDtstmTrans); 
  }
  else if (Rf_inherits(R_params, "params_mlogit_list")){
    mod = new mlogit_list(R_CohortDtstmTrans); 
  }
  else{
    Rcpp::stop("The selected statistical model is not available.");
  }
  std::unique_ptr<trans_model> uptr(mod);
  return uptr;
}

/***************************************************************************//** 
 * @ingroup dtstm
 * Simulate a Markov chain.
 * @param trans_model The transition model, which is a pointer to trans_model.  
 * @param sample A randomly sampled parameter set.
 * @param times Times at which to simulate the Markov chain.
 * @return A matrix with rows indexing @c times and columns indexing health 
 * states.
 ******************************************************************************/ 
inline arma::mat sim_markov_chain(trans_model *mod,
                                  const int &sample,
                                  const std::vector<double> &times){
  int n_states = mod->get_n_states();
  int n_times = times.size();
  arma::mat x(n_times, n_states);
  x.row(0) = mod->start_stateprobs_;
  for (int t = 1; t < n_times; ++t){
    x.row(t) = x.row(t - 1) * mod->tpmatrix(sample, times[t]);
  } 
  return x;
}

} // end namespace dtstm

} // end namespace hesim

# endif
