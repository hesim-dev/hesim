// [[Rcpp::interfaces(r, cpp)]]
#include "survival.h"
using namespace Rcpp;

/********************
* Simulating survival 
********************/
Survival::Survival(std::string dist_name, Rcpp::List surv_coefs, 
                   Rcpp::List surv_X, double r_health, double r_costs){
  dist_name_ = dist_name;
  surv_coefs_ = convert_distribution_parameters(dist_name, surv_coefs);
  surv_X_ = list_to_vecmats(surv_X);
  r_health_ = r_health;
  r_costs_ = r_costs;
  nsims_ = surv_coefs_[0].n_rows;
  nindivs_ = surv_X_[0].n_rows;
  npars_ = surv_coefs_.size();
  index_ = 0;
}

void Survival::set_index(int index){
  index_ = index;
}

std::vector<double> Survival::predict_surv_pars(int sim, int id){
  std::vector<double> y(npars_);
  for (int j = 0; j < npars_; ++j){
    y[j] = arma::dot(surv_coefs_[j].row(sim), surv_X_[j].row(id));
  }
  return y;
}

void Survival::summary1(int sim, int id, std::string type, Distribution * dist, 
                        arma::mat &output, std::vector<double> xvec){
  int J = xvec.size();
  for (int j = 0; j < J; ++j){
    output(index_, 0) = sim;
    output(index_, 1) = id;
    output(index_, 2) = xvec[j];
    if (type == "quantiles"){
      output(index_, 3) = dist->quantile(xvec[j]);
    }
    else if (type == "survival"){
      output(index_, 3) = 1 - dist->cdf(xvec[j]);
    }
    else if (type == "cumhazard"){
      output(index_, 3) = dist->cumhazard(xvec[j]);
    }
    else if (type == "hazard"){
      output(index_, 3) = dist->hazard(xvec[j]);
    }
    else if (type == "rmst"){
      output(index_, 3) = rmst(dist, xvec[j], r_health_);
    }
    else{
      Rcpp::stop("Selected type is not available.");
    }
    set_index(++index_);
  }
}

void Survival::main_loop(std::string type, std::vector<double> xvec, arma::mat &output){
  for (int s = 0; s < nsims_; ++s){
    for (int i = 0; i < nindivs_; ++i){
      Distribution *dist = select_distribution(dist_name_, predict_surv_pars(s, i));
      summary1(s, i, type, dist, output, xvec);
      delete dist;
    }
  }
}

arma::mat Survival::summary(std::vector<double> xvec, std::string type){
  arma::mat output(nsims_ * nindivs_ * xvec.size(), 4);
  main_loop(type, xvec, output);
  return output;
}

// [[Rcpp::export]]
arma::mat C_survival_summary(std::string dist_name, Rcpp::List surv_coefs, 
                             Rcpp::List surv_X, double r,
                             std::vector<double> x, std::string type){
  Survival surv(dist_name, surv_coefs, surv_X, r, 0.0);
  return surv.summary(x, type);
}


