#include <hesim/statmods/obs_index.h>

/**
 * \ingroup test
 */
class test_time_fun : public hesim::time_fun {
private:
public:
  int agecol_;
  test_time_fun(Rcpp::List L){
    agecol_ = L["agecol"];
   }
   void eval(arma::mat &m, int i, int time_old, int time_new, int state) const {
     m(i, agecol_) = m(i, agecol_) + 1;
   }
};

/**
 * \ingroup test
 */
// [[Rcpp::export]]
SEXP C_test_xptr_test_time_fun(Rcpp::List L) {
  return Rcpp::XPtr<hesim::time_fun> (new test_time_fun(L));
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
int C_test_obs_index(Rcpp::List R_input_data,
                 int strategy_id,
                 int patient_id,
                 int line = -1,
                 int health_id = -1){
  hesim::statmods::obs_index obs_index(R_input_data);
  if (line == -1 && health_id == -1){
    obs_index.set_strategy_id(strategy_id);
    obs_index.set_patient_id(patient_id);
    return obs_index();  
  }
  else if (line == -1){
    obs_index.set_strategy_id(strategy_id);
    obs_index.set_patient_id(patient_id);
    obs_index.set_health_id(health_id);
    return obs_index();  
  }
  else if (health_id == -1){
    obs_index.set_strategy_id(strategy_id);
    obs_index.set_patient_id(patient_id);
    obs_index.set_line(line);
    return obs_index();
  }
  else{
   return obs_index(strategy_id, line, patient_id, health_id);
  }
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
std::vector<int> C_test_obs_ids(Rcpp::List R_input_data, std::string member){
  hesim::statmods::obs_ids obs_ids(R_input_data);
  if (member == "strategy_id"){
    return obs_ids.strategy_id_;
  }
  else if (member == "line"){
    return obs_ids.line_;
  }
  else if (member == "patient_id"){
    return obs_ids.patient_id_;
  } 
  else if (member == "state_id"){
    return obs_ids.state_id_;
  }
  else {
    return obs_ids.transition_id_;
  }
}