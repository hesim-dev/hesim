# ifndef INPUTDATA_H
# define INPUTDATA_H
#include <RcppArmadillo.h>
#include <hesim/TimeFun.h>
#include <hesim/utils.h>

/**************
* Time function
**************/
hesim::TimeFun* get_time_fun(Rcpp::Environment R_InputData);

/**********************
* InputData base class
**********************/
class InputData {
protected:
  int strategy_index_; // Strategy index
  int mat_index_; // Matrix index including line, state, and patient
public:
  InputData(Rcpp::Environment R_InputData);
  int n_strategies_;
  std::vector<int> n_lines_;
  int n_states_;
  int n_patients_;
  int n_obs_; // Number of observations inclusive of strategies, lines, states, and patients
  int calc_n_obs();
  void set_strategy_index(int &strategy);
  void set_mat_index(int &n_lines_j, int &line, int &state, int &patient);
  void set_indices(int &strategy, int &line, int &state, int &patient);
  hesim::TimeFun* timefun_;
};

/***************************************
* InputData for single parameter models
***************************************/
class InputDataSP : public InputData{
public:
  InputDataSP(Rcpp::Environment R_InputDataSP);
  vecmats X_;
  arma::rowvec operator()() const;
  void update_X(int time_old, int time_new, int state);
};

/***************************************
* InputData for multi parameter models
***************************************/
class InputDataMP : public InputData{
public:
  InputDataMP(Rcpp::Environment R_InputDataMP);
  vecmats_2d X_;
  int param_index_;
  arma::rowvec operator()() const;
  void set_param_index(int &param);
  void update_X(int time_old, int time_new, int state);
};

// /*****************
// * InputData class
// *****************/
// class InputData {
// private:
//   int model_index_; // Model index
//   int param_index_; // Parameter index
//   int strategy_index_; // Strategy index
//   int mat_index_; // Matrix index including line, state, and patient
//   vecmats_2d init_X(SEXP X, std::string model_name);
// public:
//   InputData(Rcpp::Environment R_InputData, std::string model_name);
//   vecmats_3d X_;
//   int n_strategies_;
//   std::vector<int> n_lines_;
//   int n_states_;
//   int n_patients_;
//   int n_obs_; // Number of observations inclusive of strategies, lines, states, and patients
//   hesim::TimeFun* timefun_;
//   vecmats_3d split(vecmats_2d V);
//   void set_model_index(int &model);
//   void set_param_index(int &param);
//   void set_strategy_index(int &strategy);
//   void set_mat_index(int &n_lines_j, int &line, int &state, int &patient);
//   void set_indices(int &strategy, int &line, int &state, int &patient);
//   arma::rowvec operator()() const;
//   void update_X(int time_old, int time_new, int state);
// };

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


