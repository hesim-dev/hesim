#include <hesim/statevals.h>
#include <hesim/dtstm.h>

/***************************************************************************//** 
 * @ingroup dtstm
 * Normalized transition probabilities.
 * Normalize transition probabilities so that the probabilities of each row
 * sum to 1
 * @@param value the (possibly) unnormalized transition probabilities.
 ******************************************************************************/ 
// [[Rcpp::export]]
void C_normalize_transprobs(arma::cube &value){
  int n_slices = value.n_slices;
  int n_states = value.slice(0).n_rows;
  for (int i = 0; i < n_slices; ++i){
    for (int j = 0; j < n_states; ++j){
      double sum_probs = arma::sum(value.slice(i).row(j));
      value.slice(i).row(j) = value.slice(i).row(j)/sum_probs;
    }
  }
}

/***************************************************************************//** 
 * @ingroup dtstm
 * Simulate health state probabilities with a cohort discrete time state 
 * transition model
 * This function is exported to @c R and used in @c CohortDtstm$sim_stateprobs().
 * @param R_CohortDtstmTrans An R object of class @c CohortDtstmTrans.
 * @param times Times at which to compute state probabilities
 * @return An R data frame of the same format as dtstm::disease_prog.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::DataFrame C_cohort_dtstm_sim_stateprobs(Rcpp::Environment R_CohortDtstmTrans,
                                               std::vector<double> times){
  // Initialize 
  std::unique_ptr<hesim::dtstm::trans_model> trans_model = hesim::dtstm::trans_model::create(R_CohortDtstmTrans);
  int n_samples = trans_model->get_n_samples();
  int n_states = trans_model->get_n_states();
  int n_times = times.size();
  int N = n_samples * 
          trans_model->obs_index_.n_strategies_ *  
          trans_model->obs_index_.n_patients_ *
          n_states*
          n_times;
  hesim::stateprobs_out out(N);

  // Compute state probabilities
  int counter = 0;
  for (int s = 0; s < n_samples; ++s){
    for (int k = 0; k < trans_model->obs_index_.n_strategies_; ++k){
      trans_model->obs_index_.set_strategy_index(k);
      for (int i = 0; i< trans_model->obs_index_.n_patients_; ++ i){
        trans_model->obs_index_.set_patient_index(i);
        trans_model->obs_index_.set_time_index(0);
        
        // Simulate Markov chain
        arma::mat probs = hesim::dtstm::sim_markov_chain(trans_model.get(), s, times);
        
        // Store output
        for (int h = 0; h < n_states; ++h){
          for (int t = 0; t < n_times; ++t){
            out.sample_[counter] = s;
            out.strategy_id_[counter] = trans_model->obs_index_.get_strategy_id();
            out.patient_id_[counter] = trans_model->obs_index_.get_patient_id();
            out.grp_id_[counter] = trans_model->obs_index_.get_grp_id();
            out.patient_wt_[counter] = trans_model->obs_index_.get_patient_wt();
            out.state_id_[counter] = h;
            out.t_[counter] = times[t];
            out.prob_[counter] = probs(t, h);
            ++counter;                   
          } // end cycles loop
        } // end state loop
      } // end patient loop
    } // end strategy loop
  } // end samples loop
  
  return(out.create_R_data_frame());    
}





