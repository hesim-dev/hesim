#include <hesim/ctstm/ctstm.h>

// [[Rcpp::export]]
int tmax_max(arma::mat m){
  return m.max();
}

// [[Rcpp::export]]
Rcpp::DataFrame C_ctstm_summary(Rcpp::Environment R_CtstmTrans,
                                     std::vector<double> t,
                                     std::string type) {
  std::unique_ptr<hesim::ctstm::transmod> transmod = hesim::ctstm::transmod::create(R_CtstmTrans);
  int n_samples = transmod->get_n_samples();
  int n_strategies = transmod->obs_index_.n_strategies_;
  int n_patients = transmod->obs_index_.n_patients_;
  int n_transitions = transmod->get_n_transitions();
  int N = n_strategies * n_transitions * n_samples * n_patients * t.size();
  hesim::ctstm::transmod_summary out(N);
  
   // Loop
  int counter = 0;
  for (int s = 0; s < n_samples; ++s) { // samples loop
    for (int k = 0; k < n_strategies; ++k) { // strategies loop
      transmod->obs_index_.set_strategy_index(k);
      for (int i = 0; i < n_patients; ++i) { // patients loop
        transmod->obs_index_.set_patient_index(i);
        for (int h = 0; h < n_transitions; ++h) { // transitions loop
          
          // summarize multi-state model
          std::vector<double> res_vec;
          if (type == "cumhazard"){
            res_vec = transmod->summary(h, s, t, "cumhazard");
          }
          else if (type == "hazard"){
            res_vec = transmod->summary(h, s, t, "hazard");
          }
          
          // store results
          for (int l = 0; l < (int) t.size(); ++l){
            out.transition_id_[counter] = h;
            out.sample_[counter] = s;
            out.strategy_id_[counter] = transmod->obs_index_.get_strategy_id();
            out.patient_id_[counter] = transmod->obs_index_.get_patient_id();
            out.grp_id_[counter] = transmod->obs_index_.get_grp_id();
            out.t_[counter] = t[l];
            out.value_[counter] = res_vec[l];
            ++counter;
          }
            
        } // end transitions loop
      } // end patient loop
    } // end strategies loop
        
  } // end samples loop

   // Return
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["transition_id"] = out.transition_id_,
    Rcpp::_["sample"] = out.sample_,
    Rcpp::_["strategy_id"] = out.strategy_id_,
    Rcpp::_["patient_id"] = out.patient_id_,
    Rcpp::_["grp_id"] = out.grp_id_,
    Rcpp::_["t"] = out.t_,
    Rcpp::_["value"] = out.value_,
    Rcpp::_["stringsAsFactors"] = false
  );
  return out_df;
}
