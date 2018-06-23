# ifndef PARAMS_H
# define PARAMS_H
#include <RcppArmadillo.h>
#include <hesim/utils.h>

/*************
* Linear model
*************/
class ParamsLm   {
public:
  int sample_;
  int n_samples_;
  arma::mat coefs_;
  std::vector<double> sigma_;
  ParamsLm(Rcpp::List R_ParamsLm);
  arma::rowvec get_coefs() const;
};

/****************
* Survival models
****************/
struct SplineSurvAux{
  std::vector<double> knots_;
  std::string scale_;
  std::string timescale_;
  SplineSurvAux(Rcpp::List R_ParamsSurv);
};

struct FracPolyAux{
  std::vector<double> powers_;
  FracPolyAux(Rcpp::List R_ParamsSurv);
};

// Single survival model
class ParamsSurv   {
public:
  int sample_;
  int n_samples_;
  int n_pars_;
  vecmats coefs_;
  std::string dist_name_;
  SplineSurvAux spline_aux_;
  FracPolyAux fracpoly_aux_;
  ParamsSurv(Rcpp::List R_ParamsSurv);
  ParamsSurv(vecmats coefs, std::string dist_name, int n_samples,
             SplineSurvAux spline_aux, FracPolyAux fracpoly_aux);
};

// List of survival models
class ParamsSurvList{
public:
  std::vector<ParamsSurv> params_list_;
  int n_samples_;
  int n_models_;
  ParamsSurvList(Rcpp::List R_ParamsSurvList);
};

// Joined survival models
class ParamsJoinedSurvs {
public: 
  std::vector<ParamsSurv> params_;
  std::vector<double> times_;
  int n_samples_;
  int n_times_;
  ParamsJoinedSurvs(Rcpp::List R_ParamsJoinedSurvs);
};

# endif


