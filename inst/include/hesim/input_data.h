# ifndef INPUTDATA_H
# define INPUTDATA_H
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
 * An input data class.
 * Contains data for use as an input for prediction or sampling in a health 
 * economic simulation model. This class is instantiated using data from the R 
 * class @c hesim::input_data.
 ******************************************************************************/ 
class input_data {
private:
  int model_; ///< An integer denoting the statistical model in a list of models.
  int param_id_; ///< Parameter id used to select observation.
  int strategy_id_; ///< Strategy id used to select observation.
  int line_; ///< Line used to select observation.
  int patient_id_; ///< Patient id used to select observation.
  int health_id_; ///< Health id used to select observation.
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
  int init_n_obs_() {
    n_obs_ = 0;
    for (int i = 0; i <n_strategies_; ++i){
      n_obs_ += n_lines_[i] * n_healthvals_ * n_patients_;
    }
  }
  
  /** 
   * Initialize @c X_.
   * Initialize @c X_, a 2D vector of input matrices. There is a matrix for each
   * fitted statistical model and each parameter of a given statistical model.
   * @param[in] X The argument @c X from the @c R object @c input_data.
   * @param[out] X_ The private member @c X_.
   */    
  void init_X_(SEXP X) {
    vecmats_2d V;
    Rcpp::List L = Rcpp::as<Rcpp::List>(X);
    // If X is a list of matrices (i.e., 1 model, N parameters)
    if(Rf_isMatrix(L[0])){
      V.push_back(Rcpp::as<vecmats>(L));
    }
    // If X is a list of lists of matrices (i.e., N models, N parameters)
    else{
      V = Rcpp::as<vecmats_2d>(L);
    }
    X_ = V;
  }
public:
  /** 
   * The constructor.
   * Instantiates an input data object. 
   */  
  input_data(Rcpp::List R_input_data){
    n_strategies_ = Rcpp::as<int> (R_input_data["n_strategies"]);
    init_n_lines_(R_input_data, n_strategies_);
    init_n_healthvals_(R_input_data);
    n_patients_ = Rcpp::as<int> (R_input_data["n_patients"]);
    init_n_obs_();
    init_cum_strategy_sizes_();
    model_ = 0;
    param_id_ = 0;
    strategy_id_ = 0;
    line_ = 0;
    patient_id_ = 0;
    health_id_ = 0;
    obs_ = 0;
    timefun_ = get_time_fun(R_input_data);
    init_X_(R_input_data["X"]);
  }
  
  vecmats_2d X_; ///< Vector of input matrices. One vector for models and one for parameters.
  int obs_; ///< Observation from @c X_[model_][param_id_] based on 
            ///< @c strategy_id, @c line, @c patient_id, and @c health_id (e.g., state_id, transition_id)
  int n_strategies_; ///< Number of treatment strategies.
  std::vector<int> n_lines_; ///< Number of treatment lines for each treatment strategy.
  int n_healthvals_; ///< Number of unique health values (i.e., states, transitions).
  int n_patients_; ///< Numbber of unique patients.
  int n_obs_; ///< Number of observations inclusive of strategies, lines, health values, and patients.
  time_fun* timefun_; ///< A pointer to @c time_fun.
 
  /** 
   * Set the model
   */  
  void set_model(int model){
     model_ = model;
  }
  
  /** 
   * Set the parameter id.
   */  
  void set_param_id(int param_id) {
    param_id_ = param_id;
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
   * Set the observation number.
   * Set the row index of the matrices in @c X_ (i.e., @c obs_) based on the current values of
   * @c strategy_id_, @c line_, @c patient_id_, and @c health_id_. 
   * @param[out] obs_ The private member @c obs_.
   * @return None.
   */   
  void set_obs() {
    int strategy_row = line_ * n_patients_ * n_healthvals_ +
                     patient_id_ * n_healthvals_ +
                     health_id_;
    obs_ = strategy_row + cum_strategy_size_;
  }
  
  /** 
   * Set the observation number.
   * Set the row index of the matrices in @c X_ (i.e., @c obs_) based on new values of
   * @c strategy_id, @c line, @c patient_id, and @c health_id. Updates the private 
   * members @c strategy_id_, @c line_, @c patient_id_, and @c health_id_. 
   * @param[out] strategy_id_ The private member @c strategy_id_.
   * @param[out] line_ The private member @c line_.
   * @param[out] patient_id_ The private member @c patient_id_.
   * @param[out] health_id_ The private member @c health_id_.
   * @param[out] obs_ The private member @c obs_.
   * @return None.
   */   
  void set_obs(int strategy_id, int line, int patient_id, int health_id) {
    set_strategy_id(strategy_id);
    set_line(line);
    set_patient_id(patient_id);
    set_health_id(health_id);
    set_obs();
  }
  
  /** 
   * Return input vector.
   * The input vector associated with the current observation number (@c obs_),
   * parameter id (@c param_id_), and model (@c mdodel_).
   * @return A vector of explanatory variables.
   */  
  arma::rowvec operator()() const {
    return X_[model_][param_id_].row(obs_);
  }
};

} // end namespace statmods

} // end namespace hesim


# endif


