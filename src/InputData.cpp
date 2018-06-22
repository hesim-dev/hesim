#include <hesim/InputData.h>

// [[Rcpp::export]]
SEXP C_test_xptr_TimeFunTest(Rcpp::List L) {
  return Rcpp::XPtr<hesim::TimeFun> (new TimeFunTest(L));
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