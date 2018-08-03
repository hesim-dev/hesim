# ifndef HESIM_STATS_PARAMS_H
# define HESIM_STATS_PARAMS_H
#include <RcppArmadillo.h>
#include <hesim/utils.h>

namespace hesim {

namespace statmods {

/***************************************************************************//** 
 * Parameters of a linear model.
 ******************************************************************************/ 
class params_lm   {
public:
  int sample_;
  int n_samples_; ///< Number of random samples of the parameters.
  arma::mat coefs_; ///< A matrix of covariates where each row is a random draw
                   ///< of the parameters from the posterior distribution and
                   ///< each column is an explanatory variable.
  std::vector<double> sigma_; ///< A vector of samples of the standard error 
                             ///< of the regression model.
                           
  /** 
   * The constructor.
   * Instantiates the parameters of a linear model.
   * @param R_params_lm An object of class "params_lm" passed from the @c hesim
   * @c R package. 
   */ 
  params_lm(Rcpp::List R_params_lm) {
    coefs_ = Rcpp::as<arma::mat>(R_params_lm["coefs"]);
    sigma_ = Rcpp::as<std::vector<double> >(R_params_lm["sigma"]);
    sample_ = 0;
    n_samples_ = Rcpp::as<int> (R_params_lm["n_samples"]);
  }
};

/*******************************************************************************
* Survival models.
*******************************************************************************/
/** Internal details for hesim::statmods that should be ignored by external users.*/
namespace detail {

/***************************************************************************//** 
 * Auxiliary parameters for a spline survival model.
 * See stats::survspline for description of the member variables.
 ******************************************************************************/ 
struct survspline_aux{
  std::vector<double> knots_; 
  std::string scale_;
  std::string timescale_;
  
  /** 
   * The constructor.
   * Instantiates the auxiliary parameters of a spline survival model.
   * @param R_params_surv An object of class "params_surv" passed from the @c hesim
   * @c R package. 
   */ 
  survspline_aux(Rcpp::List R_params_surv) {
    std::string dist_name = Rcpp::as<std::string>(R_params_surv["dist"]);
    if (dist_name == "survspline"){
      Rcpp::List aux = Rcpp::as<Rcpp::List> (R_params_surv["aux"]);
      knots_ = Rcpp::as<std::vector<double> > (aux["knots"]);
      scale_ = Rcpp::as<std::string> (aux["scale"]);
      timescale_  = Rcpp::as<std::string> (aux["timescale"]); 
    }
  }
};

/***************************************************************************//** 
 * Auxiliary parameters for a fractional polynomial survival distribution.
 * See stats::fracpoly for description of the member variables.
 ******************************************************************************/ 
struct fracpoly_aux {
  std::vector<double> powers_;
  
  /** 
   * The constructor.
   * Instantiates the auxiliary parameters of a fractional polynomial model.
   * @param R_params_surv An object of class "params_surv" passed from the @c hesim
   * @c R package. 
   */ 
  fracpoly_aux(Rcpp::List R_params_surv) {
    std::string dist_name = Rcpp::as<std::string>(R_params_surv["dist"]);
    if (dist_name == "fracpoly"){
      Rcpp::List aux = Rcpp::as<Rcpp::List> (R_params_surv["aux"]);
      powers_ = Rcpp::as<std::vector<double> > (aux["powers"]);
    }  
  }
};

} // end namespace detail

/***************************************************************************//** 
 * Parameters of a survival model.
 ******************************************************************************/
class params_surv  {
public:
  int sample_;
  int n_samples_; ///< Number of random samples of the parameters.
  int n_pars_; ///< Number of parameters; for instance, a Weibull distribution
               ///< has 2 parameters (shape and scale).
  vecmats coefs_; ///< A vector of input matrices with one matrix for each
                 ///< parameter.
  std::string dist_name_; ///< A string denoting the name of the underlying 
                         ///< probability distribution.
  detail::survspline_aux spline_aux_; ///< Auxiliary parameters for a spline survival model. 
  detail::fracpoly_aux fracpoly_aux_; ///< Auxiliary parameters for a fractional polynomial survival model. 
  
  /** 
   * The constructor.
   * Instantiates the parameters of a survival model.
   * @param R_params_surv An object of class "params_surv" passed from the @c hesim
   * @c R package. 
   */ 
  params_surv(Rcpp::List R_params_surv)
    : spline_aux_(R_params_surv), fracpoly_aux_(R_params_surv) {
    coefs_ = hesim::detail::list_to_vec<vecmats, arma::mat>(R_params_surv["coefs"]);
    dist_name_ = Rcpp::as<std::string>(R_params_surv["dist"]);
    sample_ = 0;
    n_samples_ = Rcpp::as<int> (R_params_surv["n_samples"]);
    n_pars_ = coefs_.size();     
  }
};

/***************************************************************************//** 
 * Parameters of a list of survival models.
 ******************************************************************************/
class params_surv_list{
public:
  std::vector<params_surv> params_list_; ///< A vector where each element contains 
                                       ///< the parameters of a different survival model.
  int n_samples_; ///< Number of random samples of the parameters.
  int n_models_; ///< Number of survival models in the list
  
  /** 
   * The constructor.
   * Instantiates the parameters of a list/vector of survival models.
   * @param R_params_surv_list An object of class "params_surv_list" passed from the @c hesim
   * @c R package. 
   */ 
  params_surv_list(Rcpp::List R_params_surv_list) {
    n_models_ = R_params_surv_list.size();
    params_list_.reserve(n_models_);
    for (int i = 0; i < n_models_; ++i){
      Rcpp::List params_surv_i = Rcpp::as<Rcpp::List> (R_params_surv_list[i]);
      params_list_.push_back(params_surv(params_surv_i));
    }
    n_samples_ = params_list_[0].n_samples_;    
  }
};

/***************************************************************************//** 
 * Parameters of survival models joined at specified time points.
 ******************************************************************************/
class params_joined_surv {
public: 
  std::vector<params_surv> params_;
  std::vector<double> times_;
  int n_samples_;
  int n_times_;
  params_joined_surv(Rcpp::List R_params_joined_survs) {
    Rcpp::List R_params = Rcpp::as<Rcpp::List> (R_params_joined_survs["params"]);
    times_ = Rcpp::as<std::vector<double> > (R_params_joined_survs["times"]);
    n_times_ = times_.size();
    for (int i = 0; i < R_params_joined_survs.size(); ++i){
      Rcpp::List R_params_i = R_params[i];
      params_[i] = params_surv(R_params_i);
    }
    n_samples_ = params_[0].n_samples_;    
  }
};

} // end namespace statmods

} // end namespace hesim

# endif


