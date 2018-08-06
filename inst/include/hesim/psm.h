# ifndef HESIM_PSM_H
# define HESIM_PSM_H
#include <hesim/statmods/statmods.h>
#include <hesim/math/composite.h>

namespace hesim {

/** 
 * @ingroup psm
 * Partitioned survival modeling.
 */
namespace psm {

/***************************************************************************//** 
 * An abstract base class for partitioned survival models.
 * Each child is a collection of survival models used to predict the @c N-1
 * survival curves in an N-state partitioned survival model.
 ******************************************************************************/ 
class surv_mods{ 
public:
  /** 
   * The constructor.
   * Instantiates a class containing the survival models of a partitioned survival
   * model.
   */ 
  surv_mods(Rcpp::Environment R_PsmCurves);
  virtual ~surv_mods() {}
  
  /** 
   * Create a survival model object.
   * Create a unique pointer to the surv_mods abstract base class of the survival
   * model object.
   */   
  static std::unique_ptr<surv_mods> create(Rcpp::Environment R_PsmCurves);
  
  std::vector<int> strategy_id_; ///< Strategy ids.
  std::vector<int> patient_id_; ///< Patient ids.
  int virtual get_n_models() const = 0; ///< Get the number of models (i.e., number of survival curves in the
                                       /// partitioned survival model).
  int virtual get_n_samples() const = 0; ///< Get the number of times the parameters were randomly sampled.
  int virtual get_n_obs() const = 0; ///< Get the number of observations inclusive of strategies, health states, and patients.
  
  /** 
   * Summarize the survival curves from a partitioned survival model across time periods.
   * @param model A model used to summarize on of the N-1 survival curves.
   * @param sample A randomly sampled parameter set.
   * @param obs An observation number (i.e., treatment strategy and patient combination).
   * @param t A vector of times.
   * @param type "hazard" for the hazard, "cumhazard" for the cumulative hazard,
   * "survival" for survival, and "rmst" for restricted mean survival time. 
   * @param dr Discount rate. Only used to compute restricted mean survival time.
   * @return A vector of the hazard rates, cumulative hazard rates, survival probabilities,
   * or restricted mean survival times across all times in @p t.
   */ 
  std::vector<double> virtual summary(int model, int sample, int obs, std::vector<double> t, 
                                      std::string type, double dr = 0) const = 0;
  
  /** 
   * Compute the quantiles of survival curves from a partitioned survival model.
   * @param model The model (i.e., survival curve) to summarize.
   * @param sample A randomly sampled parameter set.
   * @param obs An observation number (i.e., strategy id and patient combination).
   * @param p A vector of probabilities.
   * @return A vector of quantiles evaluated at each probability in @p p.
   */ 
  std::vector<double> virtual quantile(int model, int sample, int obs, std::vector<double> p) const = 0;
}; // end of surv_mods class


/***************************************************************************//** 
 * A "list" of survival models.
 * A list of survival models where each model in the list is used to predict
 * one of the @c N-1 curves in an N-state partitioned survival model.
 ******************************************************************************/ 
class surv_list : public surv_mods  {
private:
  hesim::statmods::params_surv_list params_; ///< The parameters of a survival model.
public:
  /** 
   * The constructor.
   * Instantiates a class in which there is a separate survival model for each
   * survival curve in the partitioned survival model.
   */ 
  surv_list(Rcpp::Environment R_PsmCurves);
  vecmats_2d X_; ///< A 2-dimensional vector of matrices. The first dimension denotes
                ///< a survival model and the second dimension denotes the parameters of
                ///< a given survival model. There is an input matrix for each survival model
                ///< and parameter combination.
  int get_n_models() const;
  int get_n_samples() const;
  int get_n_obs() const;
  std::vector<double> virtual summary(int model, int sample, int obs, std::vector<double> t, 
                                      std::string type, double dr = 0) const;
  std::vector<double> virtual quantile(int model, int sample, int obs, std::vector<double> p) const;
};

/***************************************************************************//** 
 * Data container for storing summaries of survival models.
 ******************************************************************************/ 
struct surv_summary{
  /** 
   * A default constructor.
   * Instantiates a data container for a predicted survival curve.
   */ 
  surv_summary() {};
  
  /** 
   * A constructor.
   * Instantiates a data container for a predicted survival curve where all 
   * vectors in the container are initialized to a size @c n.
   */
  surv_summary(int n) {
    curve_.resize(n);
    sample_.resize(n);
    strategy_id_.resize(n);
    patient_id_.resize(n);
    x_.resize(n);
    value_.resize(n);
  }
  
  /** 
   * Create class from R object.
   * Creates the class from an @c R object of class @c Psm.
   */
  static surv_summary create_from_R(Rcpp::Environment R_Psm) {
    surv_summary out;
    Rcpp::DataFrame R_survival = Rcpp::as<Rcpp::DataFrame>(R_Psm["survival_"]);
    
    out.curve_ = Rcpp::as<std::vector<int> >(R_survival["curve"]);
    out.sample_ = Rcpp::as<std::vector<int> >(R_survival["sample"]);
    out.strategy_id_ = Rcpp::as<std::vector<int> >(R_survival["strategy_id"]);
    out.patient_id_ = Rcpp::as<std::vector<int> >(R_survival["patient_id"]);
    out.x_ = Rcpp::as<std::vector<double> >(R_survival["t"]);
    out.value_ = Rcpp::as<std::vector<double> >(R_survival["survival"]);
    out.n_states_ = Rcpp::as<int>(R_Psm["n_states"]);
    int n_curves = out.n_states_ - 1;
    out.n_sims_ = out.value_.size()/n_curves;
    
    // R to C++ indexing
    hesim::add_constant(out.curve_, -1); 
    hesim::add_constant(out.sample_, -1); 
    return out;    
  }
  
  std::vector<int> curve_; ///< The survival curve. 
  std::vector<int> sample_; ///< A randomly sampled parameter set.
  std::vector<int> strategy_id_; ///< A treatment strategy id.
  std::vector<int> patient_id_; ///< A patient id.
  std::vector<double> x_; ///< Values at which to summarize the survival curves 
                         ///< (hazard, cumulative hazard, survival, restricted mean survival time, and 
                         ///< quantiles). Either a vector of quantiles (i.e., times) or a vector of 
                        ///< probabilities (when computing quantiles).
  std::vector<double> value_; ///< The summarized value (hazard, cumulative hazard, survival proportion,
                             ///<restricted mean survival time, or quantile).
  int n_states_; ///< Number of health states.
  int n_sims_; ///< Number of simulations for each curve inclusive of treatment 
               /// strategies, patients, times, and randomly sampled parameter sets. 
               
  /** 
   * Obtain index for each vector in the struct.
   * @param sim A simulation number.
   * @param curve  A survival curve.
   * @return The index of each vector.
   */
  int index(int sim, int curve) const {
    return curve * n_sims_ + sim;
  }
};

/***************************************************************************//** 
 * Simulate health state probabilities with a partitioned survival model.
 * @param R_psm An R object of class @c Psm
 * @return An R list with two elements: 
 *  - @c stateprobs: An @c R data frame with columns equivalent to the data members in 
 *  hesim::stateprobs_out.
 *  - @c n_crossings: An integer denoting the number of times the survival 
 *  curves cross. Equivalent to hesim::psm::stateprobs::n_crossings_.
 ******************************************************************************/ 
class stateprobs{
private:
  surv_summary survcurves_; ///< Predicted survival curves.
  int n_crossings_; ///< Number of simulations (times, parameter samples, and patients) 
                   ///< at which the survival curves cross.
                   
  /** 
   * Simulate health state probabilities for a single simulation 
   * (time, parameter sample,and patient) and health state from a vector of survival 
   * curves.
   * @param sim The simulated survival curve (i.e., element in the vector of survival curves).
   * @param state The health state.
   * @return The probability that a given patient is in a specific health state at a given point in time
   * for a given randomly sampled parameter set. 
   */                      
  double sim1(int sim, int state); 
public:
  /** 
   * The constructor.
   * Instantiates a class for computing health state probabilities from (partitioned) 
   * survival curves.
   */ 
  stateprobs(Rcpp::Environment R_psm);
  
  /** 
   * Simulate the probability of being in each health states across patients, time periods,
   * and randomly sampled parameter sets. 
   * @param R_psm An @c R object of class @Psm. Must contain the data member @c Psm$survival_ simulated
   * using @c Psm$sim_survival().
   * @return An @c R list with two elements: 
   *  - @c stateprobs: An @c R data frame with columns equivalent to the data members in 
   *  hesim::stateprobs_out.
   *  - @c n_crossings: An integer denoting the number of times the survival 
   *  curves cross. Equivalent to hesim::psm::stateprobs::n_crossings_.
   */   
  Rcpp::List sim();
};


} // end psm namespace

} // end hesim namespace

# endif
