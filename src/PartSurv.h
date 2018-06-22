# ifndef PARTSURV_H
# define PARTSURV_H
#include <hesim/stats_models.h>
#include <hesim/TimeFun.h>
#include <hesim/utils.h>
#include "integrate.h"

namespace part_surv{

// Statistical modeling for predicting state values
class StateValMods{
protected:
  int n_strategies_;
  int n_patients_;
  int n_states_;
  int obs_;
public:
  StateValMods(Rcpp::Environment R_PartSurvStateVals);
  virtual ~StateValMods() {};
  std::vector<int> state_id_;
  std::vector<int> strategy_id_;
  std::vector<int> patient_id_;
  int sample_;
  void set_sample(int sample);
  void virtual set_obs(int strategy, int patient, int state) {};
  void virtual set_obs(int obs) {};
  int virtual get_n_obs() {return 0.0;}
  int virtual get_n_samples() {return 0.0;}
  double virtual predict() const {return 0.0;}
};

class LmMod : public StateValMods{
private:
  ParamsLm params_;
public:   
  LmMod(Rcpp::Environment R_PartSurvStateVals);
  arma::mat X_;
  void set_obs(int strategy, int patient, int state);
  void set_obs(int obs);
  int get_n_obs();
  int get_n_samples();
  double predict() const;
};


// Statistical modeling for predicting survival 
class SurvMods{ 
public:
  SurvMods(Rcpp::Environment R_PartSurvCurves);
  virtual ~SurvMods() {}
  std::vector<int> strategy_id_;
  std::vector<int> patient_id_;
  int n_obs_;
  int virtual get_n_models() const {return 0.0;}
  int virtual get_n_samples() const {return 0.0;}
  int virtual get_n_obs() const {return 0.0;}
  std::vector<double> virtual summary(int model, int sample, int obs, std::vector<double> t, 
                                      std::string type, double dr = 0) const {return zeros<double>(1);}
  std::vector<double> virtual quantile(int model, int sample, int obs, std::vector<double> p) const {
    return zeros<double>(1);
  }
};

class NSurvMods : public SurvMods  {
private:
  ParamsSurvList params_;
public:
  NSurvMods(Rcpp::Environment R_PartSurvCurves);
  vecmats_2d X_;
  std::vector<double> virtual summary(int model, int sample, int obs, std::vector<double> t, 
                                      std::string type, double dr = 0) const;
  std::vector<double> virtual quantile(int model, int sample, int obs, std::vector<double> p) const;
  int get_n_models() const;
  int get_n_samples() const;
  int get_n_obs() const;
};

// class NJoinedSurvs : public Survs  {
// private:
//   ParamNJoinedSurvs params_;
//   int time_id_;
//   void time_update(double t);
//   template <typename T>
//   std::vector<double> outcomes(std::vector<double> x, T fun);
// public:
//   NJoinedSurvs(Rcpp::Environment R_Curves);
//   vecmats_3d X_;
//   std::vector<double> predict_pars() const;
//   int get_n_sims() const ;
//   int get_n_models() const;
//   void set_dist();
//   std::vector<double> hazard(std::vector<double> t);
//   std::vector<double> cumhazard(std::vector<double> t);
//   std::vector<double> survival(std::vector<double> t);
//   std::vector<double> quantiles(std::vector<double> t);
//   //std::vector<double> rmst(std::vector<double> t, double dr = 0);
// };
// 

struct SurvCurvesOut{
  SurvCurvesOut();
  SurvCurvesOut(int n); // For creating an object to fill
  static SurvCurvesOut create_from_R(Rcpp::Environment R_PartSurv);
  std::vector<int> curve_;
  std::vector<int> sample_;
  std::vector<int> strategy_id_;
  std::vector<int> patient_id_;
  std::vector<double> x_;
  std::vector<double> value_;
  int n_states_;
  int n_sims_;
  int index(int sim, int curve) const;
};

class SurvCurves{
public:
  SurvCurves(Rcpp::Environment R_PartSurvCurves);
  std::unique_ptr<SurvMods> survmods_;
  Rcpp::DataFrame summary(std::vector<double> x, std::string output, double dr = 0);
};

struct StateValsOut {
  StateValsOut(int n); // For creating an object to fill
  std::vector<int> state_id_;
  std::vector<int> sample_;
  std::vector<int> strategy_id_;
  std::vector<int> patient_id_;
  std::vector<double> value_;
};

class StateVals{
public:
  StateVals(Rcpp::Environment R_PartSurvStateVals);
  std::unique_ptr<StateValMods> statevalmods_;
  Rcpp::DataFrame predict();
};

struct StateprobsOut{
  StateprobsOut();
  StateprobsOut(int n);
  static StateprobsOut create_from_R(Rcpp::DataFrame R_stateprobs);
  std::vector<int> state_id_;
  std::vector<int> sample_;
  std::vector<int> strategy_id_;
  std::vector<int> patient_id_;
  std::vector<double> t_;
  std::vector<double> prob_;
};

class Stateprobs{
private:
  SurvCurvesOut survcurves_;
  int n_crossings_;
  double sim_probs1(int sim, int state);
public:
  Stateprobs(Rcpp::Environment R_PartSurv);
  Rcpp::List sim_probs();
};

class StateValsFunc {
private:
  double r_;
  double stateval_;
  double stateprob_;
public:
  StateValsFunc(double r, double stateval, double stateprob){
    r_= r;
    stateval_ = stateval;
    stateprob_ = stateprob;
  }
  double operator()(double t){
    return exp(-r_ * t) * stateval_ * stateprob_;
  }
};

struct AucOut {
  AucOut(int n);
  std::vector<int> state_id_;
  std::vector<int> sample_;
  std::vector<int> strategy_id_;
  std::vector<int> patient_id_;
  std::vector<double> dr_;
  std::vector<std::string> type_;
  std::vector<double> value_;
};

class Auc{
private:
  std::vector<double> times_;
  int n_times_;
  int n_states_;
  StateprobsOut stateprobs_;
  std::vector<std::unique_ptr<StateValMods> > statevalmods_;
public:
  Auc(Rcpp::Environment R_PartSurv, Rcpp::DataFrame R_stateprobs, std::string type = "utility");
  std::vector<std::unique_ptr<StateValMods> > create_statevalmods_(Rcpp::Environment R_PartSurv,
                                                     std::string type);
  Rcpp::DataFrame sim(std::vector<double> dr, std::vector<std::string> type_names);
};


} // end part_surv namespace

# endif