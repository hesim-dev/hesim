# ifndef STATMODELS_H
# define STATMODELS_H
#include <hesim/distributions.h>
#include "InputData.h"
#include "Params.h"

namespace statmods{

class LinearModel {
  LinearModel(arma::rowvec x);
  double predict(arma::rowvec x, arma::rowvec beta)
};

class Surv {
public:
  Surv(arma::rowvec x);
  arma::rowvec x_;
  std::vector<double> survival(std::vector<double> t);
};

}

# endif


