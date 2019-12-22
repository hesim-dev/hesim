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
inline time_fun* get_time_fun(Rcpp::List R_object){
 if (R_object.containsElementNamed("timefun")){
    SEXP xp = R_object["timefun"];
    if (TYPEOF(xp) == EXTPTRSXP)  {
      return Rcpp::XPtr<hesim::time_fun>(xp);
    }
    else{
      Rcpp::stop("time_fun must either not be specified or be an external pointer.");
    }
  
    // To do: allow user to pass an R function
    // else {
    //   SEXP timefun_fcall = R_object["timefun"];
    //   SEXP timefun_env = R_object["timefun_env"];
    //   timefun_ = new hesim::TimeFunR(timefun_fcall, timefun_env);
    // }
  }
  else{
    return NULL;
  }  
}

/***************************************************************************//** 
 * Get ID object
 * Get the @c R object containing ID attributes from an @c R statistical model
 * class. 
 ******************************************************************************/ 
inline Rcpp::List get_id_object(Rcpp::Environment R_object){
  if(R_object.exists("input_mats") && !Rf_isNull(R_object["input_mats"])){
    return Rcpp::as<Rcpp::List>(R_object["input_mats"]); 
  }
  else{
    return Rcpp::as<Rcpp::List>(R_object["params"]);
  }
}

/*****************
 * Time intervals
 ****************/
struct time_intervals {
  int time_index_; ///< Time index used to select observation.
  std::vector<double> time_start_; ///< Vector of unique starting times.
  std::vector<double> time_stop_; ///< Vector of unique stopping times. 
  int n_times_; ///< Number of unique time intervals.
  
  /** 
   * The constructor.
   * Instantiates a @c time_intervals object. 
   */    
  time_intervals(Rcpp::List R_object){
    // Time intervals
    if (!hesim::is_null(R_object, "time_intervals")){
      Rcpp::DataFrame time_intervals = Rcpp::as<Rcpp::DataFrame>(R_object["time_intervals"]);
      time_start_ = Rcpp::as<std::vector<double> >(time_intervals["time_start"]); 
      time_stop_ = Rcpp::as<std::vector<double> >(time_intervals["time_stop"]); 
    } else{
      time_start_.push_back(0); // A single value equal to 0.
      time_stop_.push_back(INFINITY); // A single value equal to infinity.
    }
    
    // Number of time intervals
    if (!hesim::is_null(R_object, "n_times")){
      n_times_ = R_object["n_times"];
    } else{
      n_times_ = 1;
    }    
  }
};

/***************************************************************************//** 
 * Observation index
 * A class for obtaining the row index of a specfic observation contained
 * in an input matrix for prediction or sampling. This class is instantiated using 
 * data from the @c R class @c input_data from the @c hesim package.
 ******************************************************************************/ 
class obs_index {
private:
  // The current observation index
  int index_;
  
  // Other current indices
  int strategy_index_; ///< Strategy index used to select observation.
  int line_index_; ///< Line index used to select observation.
  int patient_index_; ///< Patient index used to select observation.
  int health_index_; ///< Health index used to select observation. 
  int time_index_; ///< Time index used to select observation.
  
  std::vector<int> cum_strategy_sizes_; ///< Cumulative number of observations for first @c n-1 strategies
                                        ///< where @c n is the total number of strategies. First element is
                                        ///< equal to 0. Used to get row index since that number of rows
                                        ///< varies by treatment strategy.
  int cum_strategy_size_; ///< Cumulative size of selected @c stategy_id_; that is, @c stategy_id_
                          ///< is used to select an element from @c cum_strategy_sizes_.  
  
  // Vector of IDs
  std::vector<int> strategy_id_vec_; ///< Vector of strategy IDs.
  std::vector<int> line_vec_; ///< Vector of treatment lines.
  std::vector<int> patient_id_vec_; ///< Vector of patient IDs.
  std::vector<int> health_id_vec_; ///< Vector of health IDs.
  std::vector<int> grp_id_vec_; ///<Vector of subgroup IDs.
  std::vector<double> patient_wt_vec_; ///<Vector of patient weights.
                          
  /** 
   * Initialize @c n_lines_.
   * Initialize the number of treatment lines (@c n_lines_) from @c R object @c input_data.
   * @param R_object An @c R object containing ID attributes.
   * @param n_strategies The number of strategies based on @c n_strategies_.
   * @return None.
   */  
  void init_n_lines_(Rcpp::List R_object, int n_strategies) {
    if(!hesim::is_null(R_object, "n_lines")){
      Rcpp::DataFrame n_lines_df = Rcpp::as<Rcpp::DataFrame > (R_object["n_lines"]);
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
   * @return None.
   */    
  void init_cum_strategy_sizes_() {
    cum_strategy_sizes_.reserve(n_strategies_);
    int cum_size = 0;
    cum_strategy_sizes_.push_back(cum_size);
    for (int i = 0; i < n_strategies_ - 1; ++i){
      cum_size += n_lines_.at(i) * n_patients_ * n_healthvals_ * n_times_; 
      cum_strategy_sizes_.push_back(cum_size);
    }
  }
 
  /** 
   * Initialize @c n_healthvals_
   * @param R_object An @c R object containing the element
   *                     @c n_transitions or @c n_states.
   * @return None.
   */     
  void init_n_healthvals_(Rcpp::List R_object) {
    if(!hesim::is_null(R_object, "n_states") && !hesim::is_null(R_object, "n_transitions")){
      Rcpp::stop("'n_states' and 'n_transitions' cannot both be specified.");
    }
    else if (!hesim::is_null(R_object, "n_states")){
      n_healthvals_ = R_object["n_states"];
    }
    else if (!hesim::is_null(R_object, "n_transitions")){
      n_healthvals_ = R_object["n_transitions"];
    } 
    else{
      n_healthvals_ = 1;
    }
  }
  
  /** 
   * Initialize the number of time intervals.
   * @param R_object An @c R object containing the element @c n_times.
   * parameter to initialize.
   * @return None.
   */   
  void init_n_times_(Rcpp::List R_object) {
    if (!hesim::is_null(R_object, "n_times")){
      n_times_ = R_object["n_times"];
    } else{
      n_times_ = 1;
    }
  }  
  
  /** 
   * Initialize the number of observations.
   * @param[out] n_obs_ The private member @c n_obs_.
   * @return None.
   */   
  void init_n_obs_() {
    n_obs_ = 0;
    for (int i = 0; i <n_strategies_; ++i){
      n_obs_ += n_lines_[i] * n_healthvals_ * n_patients_ * n_times_;
    }
  }
  
  /** 
   * Set the observation index 
   * Set the observation index given current values of strategy_index_,
   * line_index_, patient_index_, and health_index_.
   * @return None.
   */     
  void set_index(){
    int strategy_row = line_index_ * n_patients_ * n_healthvals_ * n_times_ +
                       patient_index_ * n_healthvals_ * n_times_ +
                       health_index_ * n_times_ +
                       time_index_;
    index_ = strategy_row + cum_strategy_size_;
  }  
  
public:
  
  // Size of each dimension
  int n_strategies_; ///< Number of treatment strategies.
  std::vector<int> n_lines_; ///< Number of treatment lines for each treatment strategy.
  int n_healthvals_; ///< Number of unique health values (i.e., states, transitions).
  int n_patients_; ///< Number of unique patients.
  int n_times_; ///< Number of unique time intervals.
  int n_obs_; ///< Number of observations inclusive of strategies, lines, patients, health values, and time intervals.
  
  // Time intervals
  std::vector<double> time_start_; ///< Vector of unique starting times.
  std::vector<double> time_stop_; ///< Vector of unique stopping times.  
  
  /** 
   * The constructor.
   * Instantiates an @c obs_index object. 
   */  
  obs_index(Rcpp::List R_object){
    // Size of each dimension
    n_strategies_ = Rcpp::as<int> (R_object["n_strategies"]);
    init_n_lines_(R_object, n_strategies_);
    n_patients_ = Rcpp::as<int> (R_object["n_patients"]);
    init_n_healthvals_(R_object);
    init_n_times_(R_object);
    init_n_obs_();
    
    init_cum_strategy_sizes_();
    
    // Current index
    strategy_index_ = 0;
    line_index_ = 0;
    patient_index_ = 0;
    health_index_ = 0;
    time_index_ = 0;
    index_ = 0;
    
    // ID vectors
    strategy_id_vec_ = Rcpp::as<std::vector<int> >(R_object["strategy_id"]);
    if (!hesim::is_null(R_object, "line")){
      line_vec_ = Rcpp::as<std::vector<int> >(R_object["line"]);
    } 
    else{
      line_vec_ = std::vector<int>(strategy_id_vec_.size(), 0);
    }
    patient_id_vec_ = Rcpp::as<std::vector<int> >(R_object["patient_id"]);
    if (!hesim::is_null(R_object, "transition_id") &&
        !hesim::is_null(R_object, "state_id")){
      Rcpp::stop("'transition_id' and 'state_id' cannot both be specified.");
    }
    if (!hesim::is_null(R_object, "transition_id")){
     health_id_vec_ = Rcpp::as<std::vector<int> >(R_object["transition_id"]); 
    }
    if (!hesim::is_null(R_object, "state_id")){
     health_id_vec_ = Rcpp::as<std::vector<int> >(R_object["state_id"]); 
    }
    
    // Time intervals
    if (!hesim::is_null(R_object, "time_intervals")){
      Rcpp::DataFrame time_intervals = Rcpp::as<Rcpp::DataFrame>(R_object["time_intervals"]);
      time_start_ = Rcpp::as<std::vector<double> >(time_intervals["time_start"]); 
      time_stop_ = Rcpp::as<std::vector<double> >(time_intervals["time_stop"]); 
    } else{
      time_start_.push_back(0); // A single value equal to 0.
      time_stop_.push_back(INFINITY); // A single value equal to infinity.
    }
    
    // Group ID and patient weights
    if (!hesim::is_null(R_object, "grp_id")){
      grp_id_vec_ = Rcpp::as<std::vector<int> >(R_object["grp_id"]); 
    }
    if (!hesim::is_null(R_object, "patient_wt")){
      patient_wt_vec_ = Rcpp::as<std::vector<double> >(R_object["patient_wt"]); 
    } else{
      patient_wt_vec_.resize(grp_id_vec_.size(), 1.0);
    }
  } // End of constructor
  
  /** 
   * Set the strategy index.
   */    
  void set_strategy_index(int strategy_index) {
    strategy_index_ = strategy_index;
    cum_strategy_size_ = cum_strategy_sizes_.at(strategy_index_);
    set_index();
  }
  
  /** 
   * Get the strategy ID given the current indices
   */    
  int get_strategy_id(){
    return strategy_id_vec_[index_];
  }
  
  /** 
   * Set the treatment line index.
   */      
  void set_line_index(int line_index) {
    line_index_ = line_index;
    set_index();
  }
  
 /** 
   * Get the treatment line given the current indices
   */    
  int get_line(){
    return line_vec_[index_];
  }
  
  /** 
   * Set the patient index.
   */   
  void set_patient_index(int patient_index) {
    patient_index_ = patient_index;
    set_index();
  }
  
  /** 
   * Get the patient ID given the current indices
   */    
  int get_patient_id(){
    return patient_id_vec_[index_];
  }  
  
  /** 
   * Get the group ID given the current indices
   */    
  int get_grp_id(){
    return grp_id_vec_[index_];
  }   
  
  /** 
   * Get the patient weight given the current indices
   */    
  int get_patient_wt(){
    return patient_wt_vec_[index_];
  }    

  /** 
   * Set the health index.
   */     
  void set_health_index(int health_index) {
    health_index_ = health_index;
    set_index();
  }
  
  /** 
   * Get the health ID given the current indices
   */    
  int get_health_id(){
    if (health_id_vec_.size() != n_obs_){
      Rcpp::stop("There is no 'health_id' in 'input_data'.");
    }
    return health_id_vec_[index_];
  }
  
  /** 
   * Set the time index.
   */     
  void set_time_index(int time_index) {
    time_index_ = time_index;
    set_index();
  }  

  /** 
   * Get the starting time given the current time index.
   */    
  double get_time_start(){
    return time_start_[time_index_];
  }    
  
  /** 
   * Get the stopping time given the current time index.
   */    
  double get_time_stop(){
    return time_stop_[time_index_];
  }     
  
  /** 
   * The observation index.
   * Computes the row index of a matrix sorted by strategy_index_, line,
   * patient_id, and line based on the current member variable.
   * @return The index.
   */ 
  int operator()() const {
    return index_;
  }  
  
  /** 
   * The observation index.
   * Computes the row index of a matrix sorted by strategy_index_, line_,
   * patient_id, health_id, and time_start based on values pass to the function; member
   * variables are updated based on these values. Updates the values of 
   * strategy_index_, line_index_, patient_index_, health_index_, time_index_, and index_. 
   * @param strategy_index The strategy index.
   * @param line_index The treatment line index.
   * @param patient_index The patient index.
   * @param health_index The health index.
   * @param time_index The time index.
   * @return The index.
   */ 
  int operator()(int strategy_index, int line_index, int patient_index, int health_index,
                 int time_index = 0) {
    strategy_index_ = strategy_index;
    line_index_ = line_index;
    patient_index_ = patient_index;
    health_index_ = health_index;
    time_index_ = time_index;
    cum_strategy_size_ = cum_strategy_sizes_.at(strategy_index_);
    set_index();
    return index_;  
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

} // end namespace statmods

} // end namespace hesim


# endif


