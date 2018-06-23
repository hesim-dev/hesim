# ifndef STATMODELS_H
# define STATMODELS_H
#include <hesim/distributions.h>
#include <hesim/Params.h>

/*******************
* Linear model (OLS)
*******************/
class StatModLm {
public:
  arma::mat X_;
  ParamsLm params_;
  StatModLm(arma::mat X, ParamsLm params);
  double predict(int sample, int obs);
  double simulate(int sample, int obs);
};

/****************
* Survival models
****************/
// Initialize distribution pointer with dummy values for parameters and actual values for
// auxillary parameters
std::unique_ptr<hesim::Distribution> init_Distribution_ptr(ParamsSurv params_surv);

class StatModSurv {
private:
  std::vector<double> predict_pars(int sample, int obs) const;

public:
  vecmats X_;
  ParamsSurv params_;
  std::unique_ptr<hesim::Distribution> dist_;
  StatModSurv(vecmats X, ParamsSurv params);
  void set_dist(int sample, int obs);

  // Return fitted survival, cummulative hazard, hazard, or restricted mean
  // survival time at a series of times
  std::vector<double> summary(int sample, int obs, std::vector<double> t, std::string type,
                              double dr = 0);

  // Return quantile of distribution
  std::vector<double> quantile(int sample, int obs, std::vector<double> p);

};

# endif


