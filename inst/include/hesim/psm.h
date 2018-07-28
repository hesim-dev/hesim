# ifndef PARTSURV_H
# define PARTSURV_H
#include <hesim/statmods.h>
#include <hesim/utils.h>
#include <hesim/integrate.h>

namespace hesim {

/** 
 * @ingroup psm
 * Partitioned survival modeling.
 */
namespace psm {

/***************************************************************************//** 
 * An abstract base class for partitioned survival models.
 * Each child is a collection of survival models used to predict the @c N-1
 * survival curves in an N-state partitioned survival model.
 ******************************************************************************/ 
class surv_mods{ 
public:
  surv_mods(Rcpp::Environment R_PartSurvCurves);
  virtual ~surv_mods() {}
  std::vector<int> strategy_id_;
  std::vector<int> patient_id_;
  int n_obs_;
  int virtual get_n_models() const = 0;
  int virtual get_n_samples() const = 0;
  int virtual get_n_obs() const = 0;
  std::vector<double> virtual summary(int model, int sample, int obs, std::vector<double> t, 
                                      std::string type, double dr = 0) const = 0;
  std::vector<double> virtual quantile(int model, int sample, int obs, std::vector<double> p) const = 0;
}; // end of surv_mods class


/***************************************************************************//** 
 * A "list" of survival models.
 * A list of survival models where each model in the list is used to predict
 * one of the @c N-1 curves in an N-state partitioned survival model.
 ******************************************************************************/ 
class surv_list : public surv_mods  {
private:
  hesim::statmods::params_surv_list params_;
public:
  surv_list(Rcpp::Environment R_PartSurvCurves);
  vecmats_2d X_;
  std::vector<double> virtual summary(int model, int sample, int obs, std::vector<double> t, 
                                      std::string type, double dr = 0) const;
  std::vector<double> virtual quantile(int model, int sample, int obs, std::vector<double> p) const;
  int get_n_models() const;
  int get_n_samples() const;
  int get_n_obs() const;
};

/***************************************************************************//** 
 * Data container for predicted survival curves.
 * Stores @c 1 of the @c N-1 survival curves predicted using surv_curves. 
 ******************************************************************************/ 
struct surv_curves_out{
  surv_curves_out();
  surv_curves_out(int n); ///< For creating an object to fill
  static surv_curves_out create_from_R(Rcpp::Environment R_PartSurv);
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

/***************************************************************************//** 
 * Data container for predicted survival curves.
 * Stores @c 1 of the @c N-1 survival curves predicted using surv_mods. 
 ******************************************************************************/ 
class surv_curves{
public:
  surv_curves(Rcpp::Environment R_PartSurvCurves);
  std::unique_ptr<surv_mods> survmods_;
  Rcpp::DataFrame summary(std::vector<double> x, std::string output, double dr = 0);
};

class stateprobs{
private:
  surv_curves_out survcurves_;
  int n_crossings_;
  double sim_probs1(int sim, int state);
public:
  stateprobs(Rcpp::Environment R_PartSurv);
  Rcpp::List sim_probs();
};


} // end psm namespace

} // end hesim namespace

# endif
