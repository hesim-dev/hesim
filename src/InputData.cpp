#include <hesim/InputData.h>

/**************
* Time function
**************/
hesim::TimeFun* get_time_fun(Rcpp::List R_InputData){
  if (R_InputData.containsElementNamed("timefun")){
    SEXP xp = R_InputData["timefun"];
    if (TYPEOF(xp) == EXTPTRSXP)  {
      return Rcpp::XPtr<hesim::TimeFun>(xp);
    }
    else{
      Rcpp::stop("TimeFun must either not be specified or be an external pointer.");
    }
  
    // To do: allow user to pass an R function
    // else {
    //   SEXP timefun_fcall = R_InputData["timefun"];
    //   SEXP timefun_env = R_InputData["timefun_env"];
    //   timefun_ = new hesim::TimeFunR(timefun_fcall, timefun_env);
    // }
  }
  else{
    return NULL;
  }
}

// [[Rcpp::export]]
SEXP C_test_xptr_TimeFunTest(Rcpp::List L) {
  return Rcpp::XPtr<hesim::TimeFun> (new TimeFunTest(L));
}

/**********
* InputData
**********/
// Class constructor
InputData::InputData(Rcpp::List R_InputData){
  n_strategies_ = Rcpp::as<int> (R_InputData["n_strategies"]);
  n_lines_ = init_n_lines(R_InputData, n_strategies_);
  n_healthvals_ = init_n_healthvals(R_InputData);
  n_patients_ = Rcpp::as<int> (R_InputData["n_patients"]);
  n_obs_ = init_n_obs();
  cum_strategy_sizes_ = init_cum_strategy_sizes();
  model_ = 0;
  param_id_ = 0;
  strategy_id_ = 0;
  line_ = 0;
  patient_id_ = 0;
  health_id_ = 0;
  obs_ = 0;
  timefun_ = get_time_fun(R_InputData);
  X_ = init_X(R_InputData["X"]);
}

std::vector<int> InputData::init_n_lines(Rcpp::List l, int n_strategies){
  if(l.containsElementNamed("n_lines")){
    Rcpp::DataFrame n_lines_df = Rcpp::as<Rcpp::DataFrame > (l["n_lines"]);
    return Rcpp::as<std::vector<int> > (n_lines_df["N"]);
  }
  else{
    std::vector<int> n_lines(n_strategies, 1);
    return n_lines;
  }
}

int InputData::init_n_healthvals(Rcpp::List l){
  if(l.containsElementNamed("n_states") && l.containsElementNamed("n_transitions")){
    Rcpp::stop("'n_states' and 'n_transitions' cannot both be specified.");
  }
  else if (l.containsElementNamed("n_states")){
    return l["n_states"];
  }
  else if (l.containsElementNamed("n_transitions")){
    return l["n_transitions"];
  } 
  else{
    return 1;
  }
}

int InputData::init_n_obs(){
  double sum = 0;
  for (int i = 0; i <n_strategies_; ++i){
    sum += n_lines_[i] * n_healthvals_ * n_patients_;
  }
  return sum;
}

std::vector<int> InputData::init_cum_strategy_sizes(){
  std::vector<int> V(n_strategies_);
  int cum_size = 0;
  for (int i = 0; i < n_strategies_; ++i){
    cum_size += n_lines_.at(i) * n_patients_ * n_healthvals_; 
    V.push_back(cum_size);
  }
  return V;
}

vecmats_2d InputData::init_X(SEXP X){
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
  return V;
}

// Setters
void InputData::set_model(int model){
  model_ = model;
}

void InputData::set_param_id(int param_id){
  param_id_ = param_id;
};

void InputData::set_strategy_id(int strategy_id){
  strategy_id_ = strategy_id;
  cum_strategy_size_ = cum_strategy_sizes_.at(strategy_id_);
}

void InputData::set_line(int line){
  line_ = line;
}

void InputData::set_patient_id(int patient_id){
  patient_id_ = patient_id;
}

void InputData::set_health_id(int health_id){
  health_id_ = health_id;
}

void InputData::set_obs(){
  int strategy_row = line_ * n_patients_ * n_healthvals_ +
                     patient_id_ * n_healthvals_ +
                     health_id_;
  obs_ = strategy_row + cum_strategy_size_;
}

void InputData::set_obs(int patient_id, int health_id){
  int strategy_row = line_ * n_patients_ * n_healthvals_ +
                     patient_id * n_healthvals_ +
                     health_id;
  obs_ = strategy_row + cum_strategy_size_;
}

void InputData::set_obs(int line, int patient_id, int health_id){
  int strategy_row = line * n_patients_ * n_healthvals_ +
                     patient_id * n_healthvals_ +
                     health_id;
  obs_ = strategy_row + cum_strategy_size_;
}

void InputData::set_obs(int strategy_id, int line, int patient_id, int health_id){
  set_strategy_id(strategy_id);
  set_obs(line, patient_id, health_id);
}

// Operator
arma::rowvec InputData::operator()() const {
  return X_[model_][param_id_].row(obs_);
}

// Unit tests
// [[Rcpp::export]]
arma::rowvec C_test_InputData(Rcpp::List R_InputData,
                                 int param_id,
                                 int strategy_id,
                                 int patient_id,
                                 int line = -1,
                                 int health_id = -1){
   InputData input_data(R_InputData);
   input_data.set_param_id(param_id);
   if (line == -1 && health_id == -1){
     input_data.set_strategy_id(strategy_id);
     input_data.set_patient_id(patient_id);
     input_data.set_obs();
   }
   else if (line == -1){
     input_data.set_strategy_id(strategy_id);
     input_data.set_patient_id(patient_id);
     input_data.set_health_id(health_id);
     input_data.set_obs();
   }
   else if (health_id == -1){
     input_data.set_strategy_id(strategy_id);
     input_data.set_patient_id(patient_id);
     input_data.set_line(line);
     input_data.set_obs();
   }
   else{
     input_data.set_obs(strategy_id, line, patient_id, health_id); 
   }
   return(input_data());
}
  