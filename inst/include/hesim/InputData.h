# ifndef INPUTDATA_H
# define INPUTDATA_H
#include <RcppArmadillo.h>
#include <hesim/TimeFun.h>
#include <hesim/utils.h>

/**************
* Time function
**************/
hesim::TimeFun* get_time_fun(Rcpp::List R_InputData);

/**********************
* InputData class
**********************/
class InputData {
private:
  int model_; // Model number in a list of models
  int param_id_; // Parameter id used to select observation from V_
  int strategy_id_; // Strategy id used to select observation from V_
  int line_; // Line used to select observation from V_
  int patient_id_; // Patient id used to select observations from V_;
  int health_id_; // Health id used to select observations from V_;
  int mat_obs_; // Observation from V_[param_id][strategy_id] based on 
                // line, patient_id, and health_id (e.g., state_id, transition_id)
  vecmats_3d V_; // One vector for parameters and one for strategies
public:
  InputData(Rcpp::List R_InputData);
  int n_strategies_;
  std::vector<int> n_lines_;
  int n_healthvals_; // Number of unique health values (i.e., states, transitions)
  int n_patients_;
  int n_obs_; // Number of observations inclusive of strategies, lines, health values, and patients
  hesim::TimeFun* timefun_;
  int calc_n_obs();
  int get_n_healthvals(Rcpp::List l);
  std::vector<int> get_n_lines(Rcpp::List l, int n_strategies);
  vecmats split(arma::mat X);
  vecmats_2d split(vecmats V);
  vecmats_3d split(vecmats_2d V);
  void set_model(int model);
  void set_param_id(int param_id);
  void set_strategy_id(int strategy_id);
  void set_line(int line);
  void set_patient_id(int patient_id);
  void set_health_id(int health_id);
  void set_obs();
  void set_obs(int patient_id, int health_id);
  void set_obs(int line, int patient_id, int health_id);
  void set_obs(int strategy_id, int line, int patient_id, int health_id);
  arma::rowvec operator()() const;
};


class TimeFunTest : public hesim::TimeFun {
private:
public:
  int agecol_;
  TimeFunTest(Rcpp::List L){
    agecol_ = L["agecol"];
   }
   void eval(arma::mat &m, int i, int time_old, int time_new, int state) const {
     m(i, agecol_) = m(i, agecol_) + 1;
   }
};


# endif


