# ifndef HESIM_STATMODS_OBS_INDEX_H
# define HESIM_STATMODS_OBS_INDEX_H
#include <RcppArmadillo.h>
#include <hesim/time_fun.h>
#include <hesim/utils.h>

namespace hesim {

namespace statmods {

/**************
* Time function
**************/
inline time_fun* get_time_fun(Rcpp::List R_input_data){
 if (R_input_data.containsElementNamed("timefun")){
    SEXP xp = R_input_data["timefun"];
    if (TYPEOF(xp) == EXTPTRSXP)  {
      return Rcpp::XPtr<hesim::time_fun>(xp);
    }
    else{
      Rcpp::stop("time_fun must either not be specified or be an external pointer.");
    }
  
    // To do: allow user to pass an R function
    // else {
    //   SEXP timefun_fcall = R_input_data["timefun"];
    //   SEXP timefun_env = R_input_data["timefun_env"];
    //   timefun_ = new hesim::TimeFunR(timefun_fcall, timefun_env);
    // }
  }
  else{
    return NULL;
  }  
}

/***************************************************************************//** 
 * Observation index
 * A class for obtaining the row index of a specfic observation contained
 * in an input matrix for prediction or sampling. This class is instantiated using 
 * data from the @c R class @c input_data from the @c hesim package.
 ******************************************************************************/ 
class obs_index {
private:
  std::vector<int> cum_strategy_sizes_; ///< Cumulative number of observations for first @c n-1 strategies
                                        ///< where @c n is the total number of strategies. First element is
                                        ///< equal to 0. Used to get row index since that number of rows
                                        ///< varies by treatment strategy.
  int cum_strategy_size_; ///< Cumulative size of selected @c stategy_id_; that is, @c stategy_id_
                          ///< is used to select an element from @c cum_strategy_sizes_.
                          
  /** 
   * Initialize @c n_lines_.
   * Initialize the number of treatment lines from @c R object @c input_data.
   * @param[in] R_input_data The @c R object @c input_data.
   * @param[in] n_strategies The number of strategies based on @c n_strategies_.
   * @param[out] n_lines_ The private member @c n_lines_.
   */  
  void init_n_lines_(Rcpp::List R_input_data, int n_strategies) {
    if(R_input_data.containsElementNamed("n_lines")){
      Rcpp::DataFrame n_lines_df = Rcpp::as<Rcpp::DataFrame > (R_input_data["n_lines"]);
      n_lines_ =  Rcpp::as<std::vector<int> > (n_lines_df["N"]);
    }
    else{
      std::vector<int> n_lines(n_strategies, 1);
      n_lines_ = n_lines;
    }    
  }
  
  /** 
   * Initialize @c cum_strategy_sizes_
   * @param[out] cum_strategy_sizes_ The private member @c cum_strategy_sizes_.
   */    
  void init_cum_strategy_sizes_() {
    cum_strategy_sizes_.reserve(n_strategies_);
    int cum_size = 0;
    cum_strategy_sizes_.push_back(cum_size);
    for (int i = 0; i < n_strategies_ - 1; ++i){
      cum_size += n_lines_.at(i) * n_patients_ * n_healthvals_; 
      cum_strategy_sizes_.push_back(cum_size);
    }
  }
 
  /** 
   * Initialize @c n_healthvals_
   * @param[out] n_healthvals_ The private member @c n_healthvals_.
   */     
  void init_n_healthvals_(Rcpp::List R_input_data) {
    if(R_input_data.containsElementNamed("n_states") && R_input_data.containsElementNamed("n_transitions")){
      Rcpp::stop("'n_states' and 'n_transitions' cannot both be specified.");
    }
    else if (R_input_data.containsElementNamed("n_states")){
      n_healthvals_ = R_input_data["n_states"];
    }
    else if (R_input_data.containsElementNamed("n_transitions")){
      n_healthvals_ = R_input_data["n_transitions"];
    } 
    else{
      n_healthvals_ = 1;
    }
  }
  
  /** 
   * Initialize the number of observations.
   * @param[out] n_obs_ The private member @c n_obs_.
   */   
  void init_n_obs_() {
    n_obs_ = 0;
    for (int i = 0; i <n_strategies_; ++i){
      n_obs_ += n_lines_[i] * n_healthvals_ * n_patients_;
    }
  }
  
public:
  int strategy_id_; ///< Strategy id used to select observation.
  int line_; ///< Line used to select observation.
  int patient_id_; ///< Patient id used to select observation.
  int health_id_; ///< Health id used to select observation.
  
  int n_strategies_; ///< Number of treatment strategies.
  std::vector<int> n_lines_; ///< Number of treatment lines for each treatment strategy.
  int n_healthvals_; ///< Number of unique health values (i.e., states, transitions).
  int n_patients_; ///< Numbber of unique patients.
  int n_obs_; ///< Number of observations inclusive of strategies, lines, health values, and patients.
  
  /** 
   * The constructor.
   * Instantiates an input data object. 
   */  
  obs_index(Rcpp::List R_input_data){
    n_strategies_ = Rcpp::as<int> (R_input_data["n_strategies"]);
    init_n_lines_(R_input_data, n_strategies_);
    init_n_healthvals_(R_input_data);
    n_patients_ = Rcpp::as<int> (R_input_data["n_patients"]);
    init_n_obs_();
    init_cum_strategy_sizes_();
    strategy_id_ = 0;
    line_ = 0;
    patient_id_ = 0;
    health_id_ = 0;
  }
  
  /** 
   * Set the strategy id.
   */    
  void set_strategy_id(int strategy_id) {
    strategy_id_ = strategy_id;
    cum_strategy_size_ = cum_strategy_sizes_.at(strategy_id_);
  }

  /** 
   * Set the treatment line.
   */      
  void set_line(int line) {
    line_ = line;
  }
  
  /** 
   * Set the patient id.
   */   
  void set_patient_id(int patient_id) {
    patient_id_ = patient_id;
  }

  /** 
   * Set the health id.
   */     
  void set_health_id(int health_id) {
    health_id_ = health_id;
  }
  
  /** 
   * The observation index.
   * Computes the row index of a matrix sorted by strategy_id, line,
   * patient_id, and line based on the current member variable.
   * @return The index.
   */ 
  int operator()() const {
    int strategy_row = line_ * n_patients_ * n_healthvals_ +
                       patient_id_ * n_healthvals_ +
                       health_id_;
    return strategy_row + cum_strategy_size_; 
  }  
  
  /** 
   * The observation index.
   * Computes the row index of a matrix sorted by strategy_id, line,
   * patient_id, and line based on values pass to the function; member
   * variables are updated based on these values.
   * @param[in] strategy_id The strategy id.
   * @param[in] line The treatment line.
   * @param[in] patient_id The patient id.
   * @param[in] health_id The health id.
   * @param[out] strategy_id_ The member variable denoting the strategy id.
   * @param[out] line_ Tee member variable denoting the treatment line.
   * @param[out] patient_id_ The member variable denoting the patient id.
   * @param[out] health_id_ The member variable denoting the health id.
   * @return The index.
   */ 
  int operator()(int strategy_id, int line, int patient_id, int health_id) {
    strategy_id_ = strategy_id;
    line_ = line;
    patient_id_ = patient_id;
    health_id_ = health_id;
    int strategy_row = line * n_patients_ * n_healthvals_ +
                       patient_id * n_healthvals_ +
                       health_id;
    cum_strategy_size_ = cum_strategy_sizes_.at(strategy_id);
    return strategy_row + cum_strategy_size_;  
  }  
};

/***************************************************************************//** 
 * Create observation index class.
 * Create the class obs_index from a @c hesim @R model.
 * @param R_model An R based model containing data and parameters. The object
 * must contain an element names "data".
 ******************************************************************************/ 
inline obs_index create_obs_index(Rcpp::Environment R_model){
  Rcpp::List R_data = Rcpp::as<Rcpp::List>(R_model["data"]);
  obs_index obs_index(R_data);
  return obs_index;
}

/***************************************************************************//** 
 * Observation identifiers.
 * Contains @c strategy_id, @c line, @c patient_id, @c state_id, 
 * and @c transition_id from the @c R @c input_data class.
 ******************************************************************************/ 
struct obs_ids {
  std::vector<int> strategy_id_;
  std::vector<int> line_;
  std::vector<int> patient_id_;
  std::vector<int> state_id_;
  std::vector<int> transition_id_;
  
  /** 
   * A constructor.
   * Instantiates the struct from an @c R model. 
   * @param R_input_data An @c R object of class @c input_data.
   */   
  obs_ids(Rcpp::List R_input_data) {
    strategy_id_ = Rcpp::as<std::vector<int> >(R_input_data["strategy_id"]);
    if (R_input_data.containsElementNamed("line")){
      line_ = Rcpp::as<std::vector<int> >(R_input_data["line"]);
    }
    patient_id_ = Rcpp::as<std::vector<int> >(R_input_data["patient_id"]);
    if (R_input_data.containsElementNamed("transition_id")){
     transition_id_ = Rcpp::as<std::vector<int> >(R_input_data["transition_id"]); 
    }
    if (R_input_data.containsElementNamed("state_id")){
     state_id_ = Rcpp::as<std::vector<int> >(R_input_data["state_id"]); 
    }
  };
  
};

} // end namespace statmods

} // end namespace hesim


# endif


