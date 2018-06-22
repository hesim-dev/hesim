# ifndef INPUTDATA_H
# define INPUTDATA_H
#include <RcppArmadillo.h>
#include <hesim/TimeFun.h>
#include <hesim/utils.h>

/**************
* Time function
**************/
inline hesim::TimeFun* get_time_fun(Rcpp::List R_InputData){
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

/**********************
* InputData class
**********************/
class InputData {
private:
  int param_id_; // Parameter id used to select observation from V_
  int strategy_id_; // Strategy id used to select observation from V_
  int line_; // Line used to select observation from V_
  int patient_id_; // Patient id used to select observations from V_;
  int health_id_; // Health id used to select observations from V_;
  int mat_obs_; // Observation from V_[param_id][strategy_id] based on 
                // line, patient_id, and health_id (e.g., state_id, transition_id)
  vecmats_2d V_; // One vector for parameters and one for strategies
public:
  InputData(Rcpp::List R_InputData){
    n_strategies_ = Rcpp::as<int> (R_InputData["n_strategies"]);
    n_lines_ = get_n_lines(R_InputData, n_strategies_);
    n_healthvals_ = get_n_healthvals(R_InputData);
    n_patients_ = Rcpp::as<int> (R_InputData["n_patients"]);
    n_obs_ = calc_n_obs();
    param_id_ = 0;
    strategy_id_ = 0;
    line_ = 0;
    patient_id_ = 0;
    health_id_ = 0;
    mat_obs_ = 0;
    timefun_ = get_time_fun(R_InputData);
    vecmats V;
    if (Rf_isMatrix(R_InputData["X"])){
      V.push_back(Rcpp::as<arma::mat> (R_InputData["X"]));
    }
    else{
      V = Rcpp::as<vecmats>(R_InputData["X"]);  
    }
    V_ = split(V);
  }
  
  int n_strategies_;
  std::vector<int> n_lines_;
  int n_healthvals_; // Number of unique health values (i.e., states, transitions)
  int n_patients_;
  int n_obs_; // Number of observations inclusive of strategies, lines, health values, and patients
  int calc_n_obs(){
    double sum = 0;
    for (int i = 0; i <n_strategies_; ++i){
      sum += n_lines_[i] * n_healthvals_ * n_patients_;
    }
    return sum;
  }
  
  int get_n_healthvals(Rcpp::List l){
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
  
  std::vector<int> get_n_lines(Rcpp::List l, int n_strategies){
    if(l.containsElementNamed("n_lines")){
      Rcpp::DataFrame n_lines_df = Rcpp::as<Rcpp::DataFrame > (l["n_lines"]);
      return Rcpp::as<std::vector<int> > (n_lines_df["N"]);
    }
    else{
      std::vector<int> n_lines(n_strategies, 1);
      return n_lines;
    }
  }

  vecmats split(arma::mat X){
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
  
  vecmats_2d split(vecmats V){
    int n_params = V.size();
    vecmats_2d V_2d(n_params);
    for (int i = 0; i < n_params; ++i){
      V_2d[i] = split(V[i]);
    }
    return V_2d;
  }
  
  void set_param_id(int param_id){
    param_id_ = param_id;
  };
  
  void set_strategy_id(int strategy_id){
    strategy_id_ = strategy_id;
  }
  
  void set_line(int line){
    line_ = line;
  }
  
  void set_patient_id(int patient_id){
    patient_id_ = patient_id;
  }
  
  void set_health_id(int health_id){
    health_id_ = health_id;
  }
  
  void set_obs(){
    mat_obs_ = line_ * n_patients_ * n_healthvals_ +
               patient_id_ * n_healthvals_ +
               health_id_;
  }
  
  void set_obs(int patient_id, int health_id){
    mat_obs_ = line_ * n_patients_ * n_healthvals_ +
                patient_id * n_healthvals_ +
                health_id;
  }
  
  void set_obs(int line, int patient_id, int health_id){
    mat_obs_ = line * n_patients_ * n_healthvals_ +
                patient_id * n_healthvals_ +
                health_id;
  }
  
  void set_obs(int strategy_id, int line, int patient_id, int health_id){
    set_strategy_id(strategy_id);
    set_obs(line, patient_id, health_id);
  }
  
  arma::rowvec operator()() const {
    return V_[param_id_][strategy_id_].row(mat_obs_);
  }
  hesim::TimeFun* timefun_;
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


