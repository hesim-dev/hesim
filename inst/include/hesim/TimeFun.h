// Based on evaluate.h from the RcppDE package
// Allows users to either pass custom R functions or C++ functors 
// Arguments are time and health state

#ifndef hesim_TimeFun_h_
#define hesim_TimeFun_h_

#include <RcppArmadillo.h>

namespace hesim{

class TimeFun {
public:
  virtual ~TimeFun() {};
  virtual void eval(arma::mat &m, int i, int time_old, int time_new, int state) const {};
  virtual void eval(arma::rowvec &v, SEXP time) const {};
};

// class TimeFunChild : public TimeFun {
// private:
// public:
//   int agecol_;
//   TimeFunChild(Rcpp::List L){
//     agecol_ = L["agecol"];
//   }
//   void eval(arma::rowvec &v, int time) const {
//     v(agecol_) = v(agecol_) + 1;
//   }
// };

class TimeFunR : public TimeFun {
public:
  TimeFunR(SEXP fcall_, SEXP env_) : fcall(fcall_), env(env_) {}
  void eval(arma::rowvec &v, SEXP time) const {}; 
  // void eval(SEXP time) {
  //     return defaultfun(time);
  // }
private:
  SEXP fcall, env;
  double defaultfun(SEXP time) {                       // essentialy same as the old evaluate
    SEXP fn = ::Rf_lang3(fcall, time, R_DotsSymbol); // this could be done with Rcpp 
    SEXP sexp_fvec = ::Rf_eval(fn, env);            // but is still a lot slower right now
    double f_result = REAL(sexp_fvec)[0];
    if (ISNAN(f_result))
        ::Rf_error("NaN value of objective function! \nPerhaps adjust the bounds.");
    return(f_result);
  }
};

} // namespace hesim
#endif