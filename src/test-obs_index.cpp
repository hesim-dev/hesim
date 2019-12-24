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
                 int strategy_index,
                 int patient_index,
                 int health_index = -1){
  hesim::statmods::obs_index obs_index(R_input_data);
  if (health_index == -1){
    obs_index.set_strategy_index(strategy_index);
    obs_index.set_patient_index(patient_index);
    return obs_index();  
  }
  else{
   return obs_index(strategy_index, patient_index, health_index);
  }
}

/**
 * \ingroup test
 */
// [[Rcpp::export]]
int C_test_obs_index_ids(Rcpp::List R_input_data,
                                      int strategy_index,
                                      int patient_index,
                                      int health_index,
                                      std::string member){
  hesim::statmods::obs_index obs_index(R_input_data);
  obs_index.set_strategy_index(strategy_index);
  obs_index.set_patient_index(patient_index);
  obs_index.set_health_index(health_index);
  
  if (member == "strategy_id"){
    return obs_index.get_strategy_id();
  }
  else if (member == "patient_id"){
    return obs_index.get_patient_id();
  } 
  else {
    return obs_index.get_health_id();
  }
}



