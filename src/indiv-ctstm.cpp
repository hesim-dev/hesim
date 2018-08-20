#include <hesim/ctstm/indiv-ctstm.h>
#include <hesim/statevals.h>

/***************************************************************************//** 
 * @ingroup ctstm
 * Simulate disease progression (i.e., a path through a multi-state model).
 * This function is exported to @c R and used in @c CtstmTrans$sim_stateprobs() and
 * in @c IndivCtstm$sim_disase().
 * @param R_CtstmTrans An R object of class @c CtstmTrans.
 * @param start_state The starting health state.
 * @param start_ages The starting age of each patient in the simulation.
 * @param max_t The maximum time to simulate the model for.
 * @param max_age The maximum age that a patient can live to.
 * @return An R data frame of the same format as ctstm::disease_prog.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::DataFrame C_ctstm_sim_disease(Rcpp::Environment R_CtstmTrans, 
                                    int start_state,
                                    std::vector<int> start_ages, 
                                    int death_state,
                                    double max_t, double max_age){
 
  // Initialize
  std::unique_ptr<hesim::ctstm::transmod> transmod = hesim::ctstm::transmod::create(R_CtstmTrans);
  std::vector<bool> absorbing = transmod->trans_mat_.absorbing_;
  start_state = start_state - 1; // Switch from R to C++ indexing
  int n_samples = transmod->get_n_samples();
  int n_strategies = transmod->get_n_strategies(); 
  std::vector<int> n_lines = transmod->get_n_lines();
  int n_patients = transmod->get_n_patients();
  hesim::ctstm::patient patient(transmod.get(), start_ages[0], 0, start_state, max_age, max_t, death_state); 
  hesim::ctstm::disease_prog disease_prog;
  int N = 0;
  for (int i = 0; i < n_strategies; ++i){
    N+= n_samples * n_lines[i] * n_patients;
  }
  disease_prog.reserve(N);
  
   // Loop
  for (int s = 0; s < n_samples; ++s){
    for (int k = 0; k < n_strategies; ++k){
      transmod->obs_index_.set_strategy_index(k);
        for (int j = 0; j < n_lines[k]; ++j){
          transmod->obs_index_.set_line_index(j);
          for (int i = 0; i < n_patients; ++i){
            transmod->obs_index_.set_patient_index(i);
            patient.age_ = start_ages[i];
            patient.time_ = 0;
            patient.state_ = start_state;
           
            while (!absorbing[patient.state_] && patient.time_ < max_t && patient.age_ < max_age){
              int from_state = patient.state_;
              double time_start = patient.time_;
              
              // Jump to new state
              patient.jump(s); // This is the key line!
              
              // Results
              if (patient.time_ < max_t || time_start == 0){ 
                disease_prog.sample_.push_back(s);
                disease_prog.strategy_id_.push_back(transmod->obs_index_.get_strategy_id());
                disease_prog.line_.push_back(transmod->obs_index_.get_line());
                disease_prog.patient_id_.push_back(transmod->obs_index_.get_patient_id());
                disease_prog.from_.push_back(from_state);
                disease_prog.to_.push_back(patient.state_);
                disease_prog.time_start_.push_back(time_start);
                disease_prog.time_stop_.push_back(patient.time_);
                if (!absorbing[patient.state_] && patient.time_ < max_t && patient.age_ < max_age){
                  disease_prog.final_.push_back(0);
                } 
                else{
                  disease_prog.final_.push_back(1);
                }
              }
              else{
                disease_prog.time_stop_.back() = patient.time_;
              }
           } // end while loop for patient
         } // end patient loop
       } // end line loop
     } // end strategy loop
   } // end parameter sampling loop
  
  // Return 
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["sample"] = disease_prog.sample_,
    Rcpp::_["strategy_id"] = disease_prog.strategy_id_,
    Rcpp::_["line"] = disease_prog.line_,
    Rcpp::_["patient_id"] = disease_prog.patient_id_,
    Rcpp::_["from"] = disease_prog.from_,
    Rcpp::_["to"] = disease_prog.to_,
    Rcpp::_["final"] = disease_prog.final_,
    Rcpp::_["time_start"] = disease_prog.time_start_,
    Rcpp::_["time_stop"] = disease_prog.time_stop_,
    Rcpp::_["stringsAsFactors"] = false
  );
  return out_df;
}

/***************************************************************************//** 
 * @ingroup ctstm
 * Simulate health state probabilities given simulated disease progression 
 * (i.e., a path through a multi-state model) from an individual-level model. 
 * Probabilities are computed by summing state occupancy over patients.
 * @param R_disease_prog An R object of simulating disease progression generated
 * using C_ctstm_sim_disease.
 * @return A vector containing health state probabilities by randomly sampled 
 * parameter set, treatment strategy, and health state. 
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::DataFrame C_ctstm_indiv_stateprobs(Rcpp::DataFrame R_disease_prog,
                                         std::vector<double> t, int n_samples,
                                         int n_strategies, std::vector<int> unique_strategy_id,
                                         std::vector<int> strategy_index,
                                         int n_states, 
                                         int n_patients,
                                         int n_lines = 1){
  hesim::ctstm::disease_prog disease_prog(R_disease_prog);
  hesim::ctstm::stateprobs_out out(n_samples * n_strategies * n_states * t.size());
  
  for(int i = 0; i < disease_prog.time_start_.size(); ++i){
    for(int j = 0; j < t.size(); ++j){
      int state;
      if (disease_prog.final_[i] == 1 && t[j] >= disease_prog.time_stop_[i]){
        state = disease_prog.to_[i];
      }
      else{
        state = disease_prog.from_[i];
      }
      int index = disease_prog.sample_[i] * n_strategies * n_states * t.size() +
                  strategy_index[i] * n_states * t.size() + 
                  state * t.size() + // need to use to when final == 1
                  j;
      if ((t[j] >= disease_prog.time_start_[i] & t[j] < disease_prog.time_stop_[i]) ||
          (t[j] >= disease_prog.time_stop_[i] & disease_prog.final_[i] == 1)){
        out.prob_[index] = out.prob_[index] + 1;
      }
    } // end loop over times  
  } // end loop over disease progression data frame 
  
  // Convert sum to proportion
  int N = n_patients * n_lines;
  for (int i = 0; i < out.prob_.size(); ++i){
    out.prob_[i] = out.prob_[i]/N;
  }
  
  // Add identifiers
  int index = 0;
  for (int s = 0; s < n_samples; ++s){
    for (int k = 0; k < n_strategies; ++k){
      for (int h = 0; h < n_states; ++h){
        for (int r = 0; r < t.size(); ++r){
          out.sample_[index] = s;
          out.strategy_id_[index] = unique_strategy_id[k];
          out.state_id_[index] = h;
          out.t_[index] = t[r];
          ++index;
        } // end time loop
      } // end state loop
    } // end strategy loop
  } // end sample loop
  
    // Return 
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["sample"] = out.sample_,
    Rcpp::_["strategy_id"] = out.strategy_id_,
    Rcpp::_["state_id"] = out.state_id_,
    Rcpp::_["t"] = out.t_,
    Rcpp::_["prob"] = out.prob_,
    Rcpp::_["stringsAsFactors"] = false
  );
  return out_df;
}

/***************************************************************************//** 
 * @ingroup ctstm
 * Simulate weighted length of stay given simulated disease progression 
 * (i.e., a path through a multi-state model) from an individual-level model. 
 * @param R_disease_prog An R object of simulating disease progression generated
 * using C_ctstm_sim_disease.
 * @param R_StateVal An R object of class @c StateVal.
 * @return A vector of weighted length of stay in each row in R_disease_prog. These
 * values are then summed by @c patient_id using @c data.table at the @c R level
 *  in the private member function @c IndivCtstm$sim_wlos. 
 ******************************************************************************/ 
// [[Rcpp::export]]
std::vector<double> C_indiv_ctstm_wlos(Rcpp::DataFrame R_disease_prog,
                                       std::vector<int> strategy_idx,
                                       std::vector<int> patient_idx,
                                       Rcpp::Environment R_StateVal,
                                       double dr, std::string type){
  hesim::ctstm::disease_prog disease_prog(R_disease_prog);
  hesim::statevals stvals(R_StateVal);
  hesim::statmods::obs_index obs_index(Rcpp::as<Rcpp::List>(R_StateVal["data"]));
  
  int N = disease_prog.sample_.size();
  std::vector<double> wlos(N);
  for (int i = 0; i < N; ++i){
    int obs = obs_index(strategy_idx[i],
                        disease_prog.line_[i],
                        patient_idx[i],
                        disease_prog.from_[i]);
    double yhat = stvals.sim(disease_prog.sample_[i], obs, type);
    wlos[i] = hesim::pv(yhat, dr,
                        disease_prog.time_start_[i],
                        disease_prog.time_stop_[i]);
  }
  return wlos;
}



