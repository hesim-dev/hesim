#include "InputData.h"

/**************
* Time function
**************/
hesim::TimeFun* get_time_fun(Rcpp::Environment R_InputData){
  SEXP xp = R_InputData["timefun"];
  if (TYPEOF(xp) == NILSXP){
    return NULL;
  }
  else if (TYPEOF(xp) == EXTPTRSXP)  {
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
};

// [[Rcpp::export]]
SEXP C_test_xptr_TimeFunTest(Rcpp::List L) {
  return Rcpp::XPtr<hesim::TimeFun> (new TimeFunTest(L));
}


/********************************
* Split InputData design matrices
********************************/
vecmats InputData_split(arma::mat X, int n_strategies, std::vector<int> n_lines,
                        int n_states, int n_patients){
  vecmats V(n_strategies);
  int counter = 0;
  // matrix for each strategy
  for (int j = 0; j < n_strategies; ++j){
    int n_row =  n_lines[j] * n_states * n_patients;
    V[j].set_size(n_row, X.n_cols);
  
    // fill matrices
    for (int l = 0; l < n_lines[j]; ++l){
      for (int s = 0; s < n_states; ++s){
        for (int p = 0; p < n_patients; ++p){
          int index = s * n_lines[j] * n_patients +
                      l * n_patients +
                      p;
          V[j].row(index) = X.row(counter);
          ++counter;
          } // end patient loop
        } // end state loop
      } // end line loop
    } // end strategy loop
  return V;
};

vecmats_2d InputData_split(vecmats V, int n_strategies, std::vector<int> n_lines,
                        int n_states, int n_patients){
  int n_params = V.size();
  vecmats_2d V_2d(n_params); 
  for (int i = 0; i < n_params; ++i){
    V_2d[i] = InputData_split(V[i], n_strategies, n_lines,
                              n_states, n_patients);
  }
  return V_2d;
}

vecmats_3d InputData_split(vecmats_2d V, int n_strategies, std::vector<int> n_lines,
                        int n_states, int n_patients){
  int n_models = V.size();
  vecmats_3d V_3d(n_models);
  for (int i = 0; i < n_models; ++i){
    V_3d[i] = InputData_split(V[i], n_strategies, n_lines,
                              n_states, n_patients);
  }
  return V_3d;
}

/*********************
* InputData base class
*********************/
int InputData::calc_n_obs(){
  double sum = 0;
  for (int i = 0; i <n_strategies_; ++i){
    sum += n_lines_[i] * n_states_ * n_patients_;
  }
  return sum;
}

InputData::InputData(Rcpp::Environment R_InputData){
  Rcpp::DataFrame n_lines_df = Rcpp::as<Rcpp::DataFrame > (R_InputData["n_lines"]);
  n_strategies_ = Rcpp::as<int> (R_InputData["n_strategies"]);
  n_lines_ = Rcpp::as<std::vector<int> > (n_lines_df["N"]);
  n_states_ = Rcpp::as<int> (R_InputData["n_states"]);
  n_patients_ = Rcpp::as<int> (R_InputData["n_patients"]);
  n_obs_ = calc_n_obs();
  strategy_index_ = 0;
  mat_index_ = 0;
  timefun_ = get_time_fun(R_InputData);
}

void InputData::set_strategy_index(int &strategy){
  strategy_index_ = strategy;
};

void InputData::set_mat_index(int &n_lines_j, int &line, int &state, int &patient){
  mat_index_ = state * n_lines_j * n_patients_ +
               line * n_patients_ +
               patient;
};

void InputData::set_indices(int &strategy, int &line, int &state, int &patient){
  set_strategy_index(strategy);
  int n_lines_j = n_lines_[strategy];
  set_mat_index(n_lines_j, line, state, patient);
}

/*************
* InputDataSP
*************/
InputDataSP::InputDataSP(Rcpp::Environment R_InputDataSP)
  : InputData(R_InputDataSP){
  arma::mat X = Rcpp::as<arma::mat>(R_InputDataSP["X"]);
  X_ = InputData_split(X, n_strategies_, n_lines_, n_states_, n_patients_);
}

arma::rowvec InputDataSP::operator()() const {
  return X_[strategy_index_].row(mat_index_);
}

void InputDataSP::update_X(int time_old, int time_new, int state){
  timefun_->eval(X_[strategy_index_], mat_index_, time_old, time_new, state);
}

/*************
* InputDataMP
*************/
InputDataMP::InputDataMP(Rcpp::Environment R_InputDataMP)
  : InputData(R_InputDataMP){
  vecmats X = Rcpp::as<vecmats>(R_InputDataMP["X"]);
  X_ = InputData_split(X, n_strategies_, n_lines_, n_states_, n_patients_);
  param_index_ = 0;
}

void InputDataMP::set_param_index(int &param){
  param_index_ = param;
};

arma::rowvec InputDataMP::operator()() const {
  return X_[param_index_][strategy_index_].row(mat_index_);
}

void InputDataMP::update_X(int time_old, int time_new, int state){
  timefun_->eval(X_[param_index_][strategy_index_], mat_index_, time_old, time_new, state);
}

// [[Rcpp::export]]
arma::rowvec C_test_InputDataMP(Rcpp::Environment R_InputData,
                                int param_index,
                                int strategy_index,
                                int line_index,
                                int state_index,
                                int patient_index){
  InputDataMP input_data(R_InputData);
  input_data.set_param_index(param_index);
  input_data.set_indices(strategy_index, line_index, state_index, patient_index);
  return(input_data());
}

// [[Rcpp::export]]
arma::rowvec C_test_TimeFunMP(Rcpp::Environment R_InputData,
                            int time_old,
                            int time_new,
                            int state){
  InputDataMP input_data(R_InputData);
  input_data.update_X(time_old, time_new, state);
  return(input_data());
}


// /***********
// * InputData
// ***********/
// vecmats_2d InputData::init_X(SEXP X, std::string model_name){
//   if(model_name == "surv_join"){
//     return  Rcpp::as<vecmats_2d>(X);
//     }
//   else{
//     vecmats V = Rcpp::as<vecmats>(X);
//     vecmats_2d V_2d;
//     V_2d.push_back(V);
//     return V_2d;
//   }
// }
// 
// InputData::InputData(Rcpp::Environment R_InputData, std::string model_name){
//   vecmats_2d V = init_X(R_InputData["X"], model_name);
//   Rcpp::DataFrame n_lines_df = Rcpp::as<Rcpp::DataFrame > (R_InputData["n_lines"]);
//   n_strategies_ = Rcpp::as<int> (R_InputData["n_strategies"]);
//   n_lines_ = Rcpp::as<std::vector<int> > (n_lines_df["N"]);
//   n_states_ = Rcpp::as<int> (R_InputData["n_states"]);
//   n_patients_ = Rcpp::as<int> (R_InputData["n_patients"]);
//   n_obs_ = V[0][0].n_rows;
//   X_ = split(V);
//   model_index_ = 0;
//   param_index_ = 0;
//   strategy_index_ = 0;
//   mat_index_ = 0;
//   timefun_ = get_time_fun(R_InputData);
// }
// 
// vecmats_3d InputData::split(vecmats_2d V){
//   int n_models = V.size();
//   vecmats_3d V_3d(n_models);
//   for (int i = 0; i < n_models; ++i){
//     int n_params = V[i].size();
//     V_3d[i].resize(n_params, vecmats(n_strategies_));
//     for (int k = 0; k < n_params; ++k){
//       int counter = 0;
//       // matrix for each strategy
//      for (int j = 0; j < n_strategies_; ++j){
//        int n_row =  n_lines_[j] * n_states_ * n_patients_;
//        V_3d[i][k][j].set_size(n_row, V[i][k].n_cols);
//        
//        // fill matrices
//        for (int l = 0; l < n_lines_[j]; ++l){
//          for (int s = 0; s < n_states_; ++s){
//            for (int p = 0; p < n_patients_; ++p){
//              int index = s * n_lines_[j] * n_patients_ +
//                          l * n_patients_ +
//                          p;
//              V_3d[i][k][j].row(index) = V[i][k].row(counter);
//              ++counter;
//              } // end patient loop
//            } // end state loop
//          } // end line loop
//        } // end strategy loop
//      } // end parameter loop
//    } // end model loop
//    return V_3d;
//  };
// 
//  void InputData::set_model_index(int &model){
//    model_index_ = model;
//  };
// 
//  void InputData::set_param_index(int &param){
//    param_index_ = param;
//  };
// 
//  void InputData::set_strategy_index(int &strategy){
//    strategy_index_ = strategy;
//  };
// 
//  void InputData::set_mat_index(int &n_lines_j, int &line, int &state, int &patient){
//    mat_index_ = state * n_lines_j * n_patients_ +
//                 line * n_patients_ +
//                 patient;
//  };
// 
// void InputData::set_indices(int &strategy, int &line, int &state, int &patient){
//    set_strategy_index(strategy);
//    int n_lines_j = n_lines_[strategy];
//    set_mat_index(n_lines_j, line, state, patient);
//  }
// 
// arma::rowvec InputData::operator()() const {
//   return X_[model_index_][param_index_][strategy_index_].row(mat_index_);
//  }
// 
// void InputData::update_X(int time_old, int time_new, int state){
//   timefun_->eval(X_[model_index_][param_index_][strategy_index_], mat_index_, time_old, time_new, state);
//  }






