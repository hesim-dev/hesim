// [[Rcpp::interfaces(r, cpp)]]
#include "survival.h"
using namespace Rcpp;

/***************
* Free functions
***************/
std::vector<double> predict_surv_pars(int sim, int id, int npars,
                                      vecmats coefs, vecmats X){
  std::vector<double> y(npars);
  for (int j = 0; j < npars; ++j){
    y[j] = arma::dot(coefs[j].row(sim), X[j].row(id));
  }
  return y;
}

/***********************
* Survival disease model
***********************/
DisModSurv::DisModSurv(Rcpp::Environment R_DisModSurv){
  dist_name_ = Rcpp::as<std::string>(R_DisModSurv["dist_name"]);
  Rcpp::List coefs_list = Rcpp::as<Rcpp::List>(R_DisModSurv["coefs"]);
  coefs_ = convert_distribution_parameters(dist_name_, coefs_list);
  Rcpp::List X_list = Rcpp::as<Rcpp::List>(R_DisModSurv["X"]);
  X_ = list_to_vecmats(X_list);
  time_length_ = Rcpp::as<double>(R_DisModSurv["time_length"]);
};

double DisModSurv::summary1(double x, std::string type, Distribution * dist, 
                                 double discount_rate, double time_length){
  if (type == "quantiles"){
    return dist->quantile(x);
  }
  else if (type == "survival"){
    return 1 - dist->cdf(x);
  }
  else if (type == "cumhazard"){
    return dist->cumhazard(x);
  }
  else if (type == "hazard"){
    return dist->hazard(x);
  }
  else if (type == "rmst"){
    return rmst(dist, x, time_length, discount_rate);
  }
  else{
    Rcpp::stop("Selected type is not available.");
  }
}

arma::mat DisModSurv::summary(std::vector<double> xvec, std::string type, double discount_rate){
  int J = xvec.size();
  int nsims = coefs_[0].n_rows;
  int nindivs = X_[0].n_rows;
  int npars = coefs_.size();
  arma::mat output(nsims * nindivs * xvec.size(), 4);
  int index = 0;
  for (int s = 0; s < nsims; ++s){
    for (int i = 0; i < nindivs; ++i){
      std::vector<double> pars = predict_surv_pars(s, i, npars, coefs_, X_);
      Distribution *dist = select_distribution(dist_name_, pars);
      for (int j = 0; j < J; ++j){
        output(index, 0) = s;
        output(index, 1) = i;
        output(index, 2) = xvec[j];
        output(index, 3) = summary1(xvec[j], type, dist, discount_rate, time_length_);
        ++index;
      }
      delete dist;
    }
  }
  return output;
}

// [[Rcpp::export]]
arma::mat C_DisModSurv_summary(Rcpp::Environment R_DisModSurv,
                             double discount_rate,
                             std::vector<double> x,
                             std::string type){
  DisModSurv surv(R_DisModSurv);
  return surv.summary(x, type, discount_rate);
}

/*************************
* Survival decision model
*************************/
DecModSurv::DecModSurv(DisModSurv dis_mod_surv)
    :dis_mod_surv_(dis_mod_surv) {
  nsims_ = dis_mod_surv_.coefs_[0].n_rows;
  nindivs_ = dis_mod_surv_.X_[0].n_rows;
  npars_ = dis_mod_surv_.coefs_.size();
}

template <typename Func>
arma::mat DecModSurv::sim_effects(int max_t, Func utility, 
                                std::vector<double> discount_rate, 
                                std::vector<int> type){
  int n_discount_rate = discount_rate.size();
  arma::mat output(nsims_ * nindivs_ * n_discount_rate, 4);
  int index = 0;
  for (int s = 0; s < nsims_; ++s){
    for (int i = 0; i < nindivs_; ++i){
      std::vector<double> pars = predict_surv_pars(s, i, npars_, dis_mod_surv_.coefs_, dis_mod_surv_.X_);
      Distribution *dist = select_distribution(dis_mod_surv_.dist_name_, pars);
      for (int r = 0; r < n_discount_rate; ++r){
          output(index, 0) = s;
          output(index, 1) = i;
          output(index, 2) = type[r];
          utility->set_indices(s, i);
          if (type[r] == 0){
            output(index, 3) = weighted_rmst(dist, max_t, discount_rate[r], 
                                             dis_mod_surv_.time_length_, utility);
          } else{
            output(index, 3) = rmst(dist, max_t, discount_rate[r], 
                                    dis_mod_surv_.time_length_);
          }
          ++index;
      }
      delete dist;
    }
  }
  return output;
}

template <typename Func>
arma::mat DecModSurv::sim_costs(std::vector<int> t, Func cost_values, int n_components, 
                        std::vector<double> discount_rate){
  int n_discount_rate = discount_rate.size();
  arma::mat output(nsims_ * nindivs_ * n_components * n_discount_rate, 4);
  int index = 0;
  for (int s = 0; s < nsims_; ++s){
    for (int i = 0; i < nindivs_; ++i){
      std::vector<double> pars = predict_surv_pars(s, i, npars_, dis_mod_surv_.coefs_, dis_mod_surv_.X_);
      Distribution *dist = select_distribution(dis_mod_surv_.dist_name_, pars);
      for (int r = 0; r < n_discount_rate; ++r){
        for (int j = 0; j < n_components; ++j){
          output(index, 0) = s;
          output(index, 1) = i;
          output(index, 2) = j;
          cost_values->set_indices(s, i, j);
          output(index, 3) = weighted_rmst(dist, t[j], discount_rate[r], 
                 dis_mod_surv_.time_length_, cost_values);
          ++index;
        }
      }
      delete dist;
    }
  }
  return output;
}

// [[Rcpp::export]]
arma::mat C_DecModSurv_effects(Rcpp::Environment R_DisModSurv,
                           int t,
                           SEXP state_values, 
                           std::vector<double> discount_rate,
                           std::vector<int> type){
  DisModSurv dis_mod_surv(R_DisModSurv);
  StateValues *utility_fun = new UtilityDefault(state_values);
  DecModSurv dec_mod_surv(dis_mod_surv);
  arma::mat sim = dec_mod_surv.sim_effects(t, utility_fun, discount_rate, type);
  delete utility_fun;
  return sim;
}

// [[Rcpp::export]]
arma::mat C_DecModSurv_costs(Rcpp::Environment R_DisModSurv, 
                           std::vector<int> t,
                           SEXP state_values, 
                           int n_components,
                           std::vector<double> discount_rate){
  DisModSurv dis_mod_surv(R_DisModSurv);
  StateValues *costs_fun = new CostsDefault(state_values);
  DecModSurv dec_mod_surv(dis_mod_surv);
  arma::mat sim = dec_mod_surv.sim_costs(t, costs_fun, n_components, discount_rate);
  delete costs_fun;
  return sim;
}

/*************************
* Partition survival model
**************************/
PartitionSurvival::PartitionSurvival(std::string dist_name, Rcpp::List pfs_coefs, Rcpp::List pfs_X,
                                     Rcpp::List os_coefs, Rcpp::List os_X,
                                     double r_health, double r_costs, double time_length) {
  dist_name_ = dist_name;
  pfs_coefs_ = convert_distribution_parameters(dist_name, pfs_coefs);
  pfs_X_ = list_to_vecmats(pfs_X);
  os_coefs_ = convert_distribution_parameters(dist_name, os_coefs);
  os_X_ = list_to_vecmats(os_X);
  r_health_ = r_health;
  r_costs_ = r_costs;
  time_length_ = time_length;
  nsims_ = pfs_coefs_[0].n_rows;
  nindivs_ = pfs_X_[0].n_rows;
  npars_ = pfs_coefs_.size();
}

std::vector<double> PartitionSurvival::predict_pars(int s, int i, int k){
  switch(k)
  {
  case 0: 
    return predict_surv_pars(s, i, npars_, pfs_coefs_, pfs_X_);
  case 1:
    return predict_surv_pars(s, i, npars_, os_coefs_, os_X_);
  default:
      Rcpp::stop("Survival summary can only be estimated for PFS and OS.");
  }
}

arma::mat PartitionSurvival::summary(std::vector<double> xvec, std::string type){
  int J = xvec.size();
  arma::mat output(nsims_ * nindivs_ * xvec.size(), 4);
  int index = 0;
  for (int s = 0; s < nsims_; ++s){
    for (int i = 0; i < nindivs_; ++i){
      for (int k = 0; k < 2; ++k){
        std::vector<double> pars = predict_pars(s, i, k);
        Distribution *dist = select_distribution(dist_name_, pars);
        for (int j = 0; j < J; ++j){
          output(index, 0) = s;
          output(index, 1) = i;
          output(index, 2) = xvec[j];
          output(index, 3) = DisModSurv::summary1(xvec[j], type, dist, r_health_, time_length_);
          ++index;
        }
        delete dist;
      }
    }
  }
  return output;
}

