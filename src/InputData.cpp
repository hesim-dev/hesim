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
InputData::InputData(Rcpp::List R_InputData){
  n_strategies_ = Rcpp::as<int> (R_InputData["n_strategies"]);
  n_lines_ = get_n_lines(R_InputData, n_strategies_);
  n_healthvals_ = get_n_healthvals(R_InputData);
  n_patients_ = Rcpp::as<int> (R_InputData["n_patients"]);
  n_obs_ = calc_n_obs();
  model_ = 0;
  param_id_ = 0;
  strategy_id_ = 0;
  line_ = 0;
  patient_id_ = 0;
  health_id_ = 0;
  mat_obs_ = 0;
  timefun_ = get_time_fun(R_InputData);
  vecmats_2d V;
  if (Rf_isMatrix(R_InputData["X"])){
    V.resize(1);
    V[0].push_back(Rcpp::as<arma::mat> (R_InputData["X"]));
  }
  else{
    V.resize(1);
    V[0] = Rcpp::as<vecmats>(R_InputData["X"]);  
  }
  V_ = split(V);
}

int InputData::calc_n_obs(){
  double sum = 0;
  for (int i = 0; i <n_strategies_; ++i){
    sum += n_lines_[i] * n_healthvals_ * n_patients_;
  }
  return sum;
}

int InputData::get_n_healthvals(Rcpp::List l){
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

std::vector<int> InputData::get_n_lines(Rcpp::List l, int n_strategies){
  if(l.containsElementNamed("n_lines")){
    Rcpp::DataFrame n_lines_df = Rcpp::as<Rcpp::DataFrame > (l["n_lines"]);
    return Rcpp::as<std::vector<int> > (n_lines_df["N"]);
  }
  else{
    std::vector<int> n_lines(n_strategies, 1);
    return n_lines;
  }
}

vecmats InputData::split(arma::mat X){
  vecmats V(n_strategies_);
  int counter = 0;
  // matrix for each strategy
  for (int k = 0; k < n_strategies_; ++k){
    int n_row =  n_lines_[k] * n_healthvals_ * n_patients_;
    V[k].set_size(n_row, X.n_cols);

    // fill matrices
    for (int l = 0; l < n_lines_[k]; ++l){
      for (int p = 0; p < n_patients_; ++p){
        for (int h = 0; h < n_healthvals_; ++h){
          int obs = l * n_patients_ * n_healthvals_ +
                      p * n_healthvals_ +
                      h;
          V[k].row(obs) = X.row(counter);
          ++counter;
          } // end health loop
        } // end patient loop
      } // end line loop
    } // end strategy loop
  return V;
};

vecmats_2d InputData::split(vecmats V){
  int n_params = V.size();
  vecmats_2d V_2d(n_params);
  for (int i = 0; i < n_params; ++i){
    V_2d[i] = split(V[i]);
  }
  return V_2d;
}

vecmats_3d InputData::split(vecmats_2d V){
  int n_models = V.size();
  vecmats_3d V_3d(n_models);
  for (int i = 0; i < n_models; ++i){
    V_3d[i] = split(V[i]);
  }
  return V_3d;
}

void InputData::set_model(int model){
  model_ = model;
}

void InputData::set_param_id(int param_id){
  param_id_ = param_id;
};

void InputData::set_strategy_id(int strategy_id){
  strategy_id_ = strategy_id;
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
  mat_obs_ = line_ * n_patients_ * n_healthvals_ +
             patient_id_ * n_healthvals_ +
             health_id_;
}

void InputData::set_obs(int patient_id, int health_id){
  mat_obs_ = line_ * n_patients_ * n_healthvals_ +
              patient_id * n_healthvals_ +
              health_id;
}

void InputData::set_obs(int line, int patient_id, int health_id){
  mat_obs_ = line * n_patients_ * n_healthvals_ +
              patient_id * n_healthvals_ +
              health_id;
}

void InputData::set_obs(int strategy_id, int line, int patient_id, int health_id){
  set_strategy_id(strategy_id);
  set_obs(line, patient_id, health_id);
}

arma::rowvec InputData::operator()() const {
  return V_[model_][param_id_][strategy_id_].row(mat_obs_);
}

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