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
    while (time > obs_index_.get_time_stop()){
      obs_index_.set_time_index(obs_index_.get_time_index() + 1);
    }
    return value_.slice(obs_index_() + s * n_obs_);
  }

}; // end class tparams_transprobs

/** 
 * A factory function.
 * Creates a discrete time cohort transition model of the class specified
 * in @c R_CohortDtstmTrans. The class must be inherited from dtstm::trans_model.
 * @param R_CohortDtstmTrans An @c R object of class @c CohortDtstmTrans.
 * @return A unique pointer to the abstract base class dtstm::trans_model. 
 */ 
inline std::unique_ptr<trans_model> trans_model::create(Rcpp::Environment R_CohortDtstmTrans) {
  Rcpp::List R_params = Rcpp::as<Rcpp::List>(R_CohortDtstmTrans["params"]);
  trans_model * mod;
  if (Rf_inherits(R_params, "tparams_transprobs")){
    mod = new tparams_transprobs(R_CohortDtstmTrans); 
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