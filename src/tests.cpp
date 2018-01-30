#include <RcppArmadillo.h>
#include "distributions.h"
using namespace Rcpp;

RCPP_MODULE(Distributions){
  class_<Distribution>("Distribution")
  .method("pdf", &Distribution::pdf)
  .method("cdf", &Distribution::cdf)
  .method("quantile", &Distribution::quantile)
  .method("hazard", &Distribution::hazard)
  .method("cumhazard", &Distribution::cumhazard)
  .method("random", &Distribution::random)
  ;

  class_<Exponential>("Exponential")
    .derives<Distribution>("Distribution")
    .constructor<double>()
    .method("pdf", &Exponential::pdf)
    .method("cdf", &Exponential::cdf)
    .method("quantile", &Exponential::quantile)
    .method("hazard", &Exponential::hazard)
    .method("cumhazard", &Exponential::cumhazard)
    .method("random", &Exponential::random)
  ;

  class_<Weibull>("Weibull")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &Weibull::pdf)
    .method("cdf", &Weibull::cdf)
    .method("quantile", &Weibull::quantile)
    .method("hazard", &Weibull::hazard)
    .method("cumhazard", &Weibull::cumhazard)
    .method("random", &Weibull::random)
  ;
  
  class_<Gamma>("Gamma")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &Gamma::pdf)
    .method("cdf", &Gamma::cdf)
    .method("quantile", &Gamma::quantile)
    .method("hazard", &Gamma::hazard)
    .method("cumhazard", &Gamma::cumhazard)
    .method("random", &Gamma::random)
  ;
  
  class_<Lognormal>("Lognormal")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &Lognormal::pdf)
    .method("cdf", &Lognormal::cdf)
    .method("quantile", &Lognormal::quantile)
    .method("hazard", &Lognormal::hazard)
    .method("cumhazard", &Lognormal::cumhazard)
    .method("random", &Lognormal::random)
  ;
  
  class_<Gompertz>("Gompertz")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &Gompertz::pdf)
    .method("cdf", &Gompertz::cdf)
    .method("quantile", &Gompertz::quantile)
    .method("hazard", &Gompertz::hazard)
    .method("cumhazard", &Gompertz::cumhazard)
    .method("random", &Gompertz::random)
  ;
  
  class_<LogLogistic>("LogLogistic")
    .derives<Distribution>("Distribution")
    .constructor<double, double>()
    .method("pdf", &LogLogistic::pdf)
    .method("cdf", &LogLogistic::cdf)
    .method("quantile", &LogLogistic::quantile)
    .method("hazard", &LogLogistic::hazard)
    .method("cumhazard", &LogLogistic::cumhazard)
    .method("random", &LogLogistic::random)
  ;
  
  class_<GeneralizedGamma>("GeneralizedGamma")
    .derives<Distribution>("Distribution")
    .constructor<double, double, double>()
    .method("pdf", &GeneralizedGamma::pdf)
    .method("cdf", &GeneralizedGamma::cdf)
    .method("quantile", &GeneralizedGamma::quantile)
    .method("hazard", &GeneralizedGamma::hazard)
    .method("cumhazard", &GeneralizedGamma::cumhazard)
    .method("random", &GeneralizedGamma::random)
  ;
  
    class_<SurvSplines>("SurvSplines")
    .derives<Distribution>("Distribution")
    .constructor<std::vector<double>, std::vector<double>, std::string, std::string>()
    .method("linear_predict", &SurvSplines::linear_predict)
    .method("linear_predict_dx", &SurvSplines::linear_predict_dx)
    .method("pdf", &SurvSplines::pdf)
    .method("cdf", &SurvSplines::cdf)
    .method("quantile", &SurvSplines::quantile)
    .method("hazard", &SurvSplines::hazard)
    .method("cumhazard", &SurvSplines::cumhazard)
    .method("random", &SurvSplines::random)
  ;
}