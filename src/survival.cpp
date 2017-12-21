// [[Rcpp::interfaces(r, cpp)]]
#include "survival.h"
using namespace Rcpp;

/********************
* Simulating survival 
********************/
Survival::Survival(std::string dist_name, Rcpp::List surv_coefs, 
                   Rcpp::List surv_X){
  dist_name_ = dist_name;
  surv_coefs_ = convert_distribution_parameters(dist_name, surv_coefs);
  surv_X_ = list_to_vecmats(surv_X);
  nsims_ = surv_coefs_[0].n_rows;
  nindivs_ = surv_X_[0].n_rows;
  npars_ = surv_coefs_.size();
}

int Survival::get_nsims(){
  return nsims_;
}

int Survival::get_nindivs(){
  return nindivs_;
}

std::vector<double> Survival::predict_surv_pars(int sim, int id){
  std::vector<double> y(npars_);
  for (int j = 0; j < npars_; ++j){
    y[j] = arma::dot(surv_coefs_[j].row(sim), surv_X_[j].row(id));
  }
  return y;
}

template <typename Func>
arma::mat Survival::main_loop(Func fun, std::vector<double> xvec){
  int J = xvec.size();
  arma::mat output(nsims_ * nindivs_ * J, 4);
  int it = 0;
  for (int s = 0; s < nsims_; ++s){
    for (int i = 0; i < nindivs_; ++i){
      Distribution *dist = select_distribution(dist_name_, predict_surv_pars(s, i));
      for (int j = 0; j < J; ++j){
        output(it, 0) = s;
        output(it, 1) = i;
        output(it, 2) = xvec[j];
        output(it, 3) = fun(dist, xvec[j]);
        ++it;
      }
      delete dist;
    }
  }
  return output;
}

arma::mat Survival::quantiles(std::vector<double> q){
  QuantileFun qf;
  return main_loop(qf, q);
}

arma::mat Survival::hazard(std::vector<double> t){
  HazardFun hf;
  return main_loop(hf, t);
}

arma::mat Survival::cumhazard(std::vector<double> t){
  CumhazardFun chf;
  return main_loop(chf, t);
}

arma::mat Survival::surv(std::vector<double> t){
  SurvFun surv;
  return main_loop(surv, t);
}

// [[Rcpp::export]]
arma::mat C_survival_summary(std::string dist_name, Rcpp::List surv_coefs, 
                              Rcpp::List surv_X, std::vector<double> x,
                              std::string type){
  Survival surv(dist_name, surv_coefs, surv_X);
  if (type == "quantiles"){
    return surv.quantiles(x);
  }
  else if (type == "hazard"){
    return surv.hazard(x);
  }
  else if (type == "cumhazard"){
    return surv.cumhazard(x);
  }
  else if (type == "survival"){
    return surv.surv(x);
  }
  else{
    Rcpp::stop("Selected type is not available");
  }
}


