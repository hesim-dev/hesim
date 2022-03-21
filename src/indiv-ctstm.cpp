#include <hesim/ctstm/indiv-ctstm.h>
#include <hesim/statevals.h>
#include <hesim/check_R_infinity.h>

/***************************************************************************//** 
 * @ingroup ctstm
 * Simulate disease progression (i.e., a path through a multi-state model).
 * This function is exported to @c R and used in @c CtstmTrans$sim_stateprobs() and
 * in @c IndivCtstm$sim_disase().
 * @param R_CtstmTrans An R object of class @c CtstmTrans.
 * @param start_state The starting health state for each patient and random sample
 * of the parameter sets.
 * @param start_age The starting age of each patient in the simulation.
 * @param start_time The starting time of the simulation.
 * @param max_t The maximum time to simulate the model for.
 * @param max_age The maximum age that a patient can live to.
 * @param progress An integer specifying the PSA iteration (i.e.,m sample) that
 * should be printed every @c progrtess PSA iterations.
 * @param max_jumps A positive integer denoting the maximum number of jumps
 * (i.e., transitions to new states) that each patient can make. Default is -1, 
 * in which case there is no limit and the duration of the simulation is
 *  determined by @c max_t and @c max_age.
 * @param sample_start, sample_stop Starting and stopping parameter sample to 
 * iterate over. Default is to start at 0 and stop at the final parameter 
 * sample.
 * @return An R data frame of the same format as ctstm::disease_prog.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::DataFrame C_ctstm_sim_disease(Rcpp::Environment R_CtstmTrans, 
                                    std::vector<int> start_state,
                                    std::vector<double> start_age, 
                                    std::vector<double> start_time,
                                    int death_state,
                                    std::string clock, 
                                    std::vector<int> reset_states,
                                    double max_t, double max_age,
                                    int progress,
                                    int max_jumps = -1,
                                    int sample_start = 0,
                                    int sample_stop = -1){
 
  // Initialize
  hesim::check_R_infinity(max_t);
  std::unique_ptr<hesim::ctstm::transmod> transmod = hesim::ctstm::transmod::create(R_CtstmTrans);
  transmod->set_max_x(max_t);
  std::vector<bool> absorbing = transmod->trans_mat_.absorbing_;
  int n_samples = transmod->get_n_samples();
  if (sample_stop < 0) {
    sample_stop = n_samples - 1;
  }
  int n_strategies = transmod->get_n_strategies(); 
  int n_patients = transmod->get_n_patients();
  hesim::ctstm::patient patient(transmod.get(), start_age[0], start_time[0],
                                start_state[0], max_age, max_t, death_state, 
                                clock, reset_states); 
  hesim::ctstm::disease_prog disease_prog;
  int N = n_samples * n_strategies * n_patients;
  disease_prog.reserve(N);
  
   // Loop
  for (int s = sample_start; s < sample_stop + 1; ++s){
    if (progress > 0){
      if ((s + 1) % progress == 0){ // R-based indexing
        Rcpp::Rcout << "sample = " << s + 1 << std::endl;  
      }
    } 
    for (int k = 0; k < n_strategies; ++k){
      transmod->obs_index_.set_strategy_index(k);
        for (int i = 0; i < n_patients; ++i){
          transmod->obs_index_.set_patient_index(i);
          patient.age_ = start_age[i];
          patient.time_ = start_time[i];
          patient.clockmix_time_ = start_time[i];
          patient.state_ = start_state[i];
          patient.max_t_ = max_t;
         
         int jump = 1; 
          while (!absorbing[patient.state_] && patient.time_ < max_t && 
                 patient.age_ < max_age &&
                 (max_jumps < 0 || jump <= max_jumps)){
            int from_state = patient.state_;
            double time_start = patient.time_;
            
            // Jump to new state
            patient.jump(s); // This is the key line!
            
            // Results
            disease_prog.sample_.push_back(s);
            disease_prog.strategy_id_.push_back(transmod->obs_index_.get_strategy_id());
            disease_prog.patient_id_.push_back(transmod->obs_index_.get_patient_id());
            disease_prog.grp_id_.push_back(transmod->obs_index_.get_grp_id());
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
            ++jump;
          
         } // end while loop for patient
       } // end patient loop
     } // end strategy loop
   } // end parameter sampling loop
  
  // Return 
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["sample"] = disease_prog.sample_,
    Rcpp::_["strategy_id"] = disease_prog.strategy_id_,
    Rcpp::_["patient_id"] = disease_prog.patient_id_,
    Rcpp::_["grp_id"] = disease_prog.grp_id_,
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
                                         int n_strategies, 
                                         std::vector<int> unique_strategy_id,
                                         std::vector<int> strategy_index,
                                         int n_grps,
                                         std::vector<int> unique_grp_id,
                                         std::vector<int> grp_index,
                                         int n_states, 
                                         int n_patients){
  hesim::ctstm::disease_prog disease_prog(R_disease_prog);
  hesim::stateprobs_out out(n_samples * n_strategies * n_grps * n_states * t.size());
  
  for(int i = 0; i < disease_prog.time_start_.size(); ++i){
    for(int j = 0; j < t.size(); ++j){
      int state;
      if (disease_prog.final_[i] == 1 && t[j] >= disease_prog.time_stop_[i]){
        state = disease_prog.to_[i];
      }
      else{
        state = disease_prog.from_[i];
      }
      int index = disease_prog.sample_[i] * n_strategies * n_grps * n_states * t.size() +
                  strategy_index[i] * n_grps * n_states * t.size() + 
                  grp_index[i] * n_states * t.size() + 
                  state * t.size() + // need to use to when final == 1
                  j;
      if ((t[j] >= disease_prog.time_start_[i] && t[j] < disease_prog.time_stop_[i]) ||
          (t[j] >= disease_prog.time_stop_[i] && disease_prog.final_[i] == 1)){
        out.prob_[index] = out.prob_[index] + 1;
      }
    } // end loop over times  
  } // end loop over disease progression data frame 
  
  // Convert sum to proportion
  for (int i = 0; i < out.prob_.size(); ++i){
    out.prob_[i] = out.prob_[i]/n_patients;
  }
  
  // Add identifiers
  int index = 0;
  for (int s = 0; s < n_samples; ++s){
    for (int k = 0; k < n_strategies; ++k){
      for (int g = 0; g < n_grps; ++g){
        for (int h = 0; h < n_states; ++h){
          for (int r = 0; r < t.size(); ++r){
            out.sample_[index] = s;
            out.strategy_id_[index] = unique_strategy_id[k];
            out.grp_id_[index] = unique_grp_id[g];
            out.state_id_[index] = h;
            out.t_[index] = t[r];
            ++index;
          } // end time loop
        } // end state loop
      } // end subgroup loop
    } // end strategy loop
  } // end sample loop
  
    // Return 
  Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
    Rcpp::_["sample"] = out.sample_,
    Rcpp::_["strategy_id"] = out.strategy_id_,
    Rcpp::_["grp_id"] = out.grp_id_,
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
 * @param strategy_idx The strategy index starting at 0.
 * @param patient_idx The patient index starting at 0.
 * @param R_StateVal An R object of class @c StateVal.
 * @param dr The discount rate
 * @param type @c predict for mean values or @c random for random samples.
 * @param max_time Maximum time duration to compute costs once a patient has
 * entered a (new) healthy state. 
 * @return A vector of weighted length of stay in each row in R_disease_prog. These
 * values are then summed by @c patient_id using @c data.table at the @c R level
 *  in the private member function @c IndivCtstm$sim_wlos. 
 ******************************************************************************/ 
// [[Rcpp::export]]
std::vector<double> C_indiv_ctstm_wlos(Rcpp::DataFrame R_disease_prog,
                                       std::vector<int> strategy_idx,
                                       std::vector<int> patient_idx,
                                       Rcpp::Environment R_StateVal,
                                       double dr, std::string type,
                                       double max_time){
  hesim::ctstm::disease_prog disease_prog(R_disease_prog);
  hesim::statevals stvals(R_StateVal);
  hesim::statmods::obs_index obs_index(hesim::statmods::get_id_object(R_StateVal));
  hesim::check_R_infinity(max_time);
  bool time_reset = Rcpp::as<bool>(R_StateVal["time_reset"]);
  
  int N = disease_prog.sample_.size();
  std::vector<double> wlos(N);
  int time_index = 0;
  double time_stop_it;

  for (int i = 0; i < N; ++i){
    double wlos_it = 0; //Weighted LOS by row and time interval
    double time_start_it = disease_prog.time_start_[i];
    double time_stop_i_max = disease_prog.time_stop_[i];
    if (!std::isinf(max_time)){
        time_stop_i_max = std::min(disease_prog.time_stop_[i],
                                   disease_prog.time_start_[i] + max_time);
    }
    if (i > 0){
      if (time_reset){   // If time_reset = TRUE, then reset time when a patient enters a new state
        time_index = 0;  // In disease_prog_ each row denotes a transition to a new health state
      }
      else{ // Otherwise, reset time for each new patient
        if (disease_prog.patient_id_[i] != disease_prog.patient_id_[i-1] ||
            disease_prog.sample_[i] != disease_prog.sample_[i-1]){
              time_index = 0;
        }
      }
    } 
    int t_start = time_index;
    for (int t = t_start; t < obs_index.n_times_; ++t){
      int obs = obs_index(strategy_idx[i],
                          patient_idx[i],
                          disease_prog.from_[i],
                          t);  
      
      // Get stopping time
      if (time_reset){
        // If time_reset == true
        if (obs_index.get_time_start() > disease_prog.time_stop_[i] - disease_prog.time_start_[i]){
          break;
        }   
        time_stop_it = std::min(disease_prog.time_start_[i] + obs_index.get_time_stop(), 
                                time_stop_i_max);
      }
      else{
        // If time_reset == false
        if (obs_index.get_time_start() > disease_prog.time_stop_[i]){
          break;
        }     
        time_stop_it = std::min(obs_index.get_time_stop(), 
                                time_stop_i_max);
      }
      
      // Compute present value
      double yhat = stvals.sim(disease_prog.sample_[i], obs, type);
      wlos_it = hesim::pv(yhat, dr,
                          time_start_it,
                          time_stop_it);
      wlos[i] += wlos_it; 
      
      // Update starting time and time index
      time_start_it = time_stop_it;
      if (t < (obs_index.n_times_ - 1) && time_stop_it >= obs_index.get_time_stop()){
        time_index = t + 1; 
      }
    } // end loop over time intervals
  } // end loop over rows
  return wlos;
}

/***************************************************************************//** 
 * @ingroup ctstm
 * One-time expected values.
 * Simulate expected values in an individual patient simulation based on one-time 
 * state values that occur when a patient enters a health state. Values 
 * accrue each time a patient enters a new state and are discounted based on time 
 * of entrance into that state. 
 * @param R_disease_prog An R object of simulating disease progression generated
 * using C_ctstm_sim_disease.
 * @param strategy_idx The strategy index starting at 0.
 * @param patient_idx The patient index starting at 0.
 * @param R_StateVal An R object of class @c StateVal.
 * @param dr The discount rate.
 * @param type @c predict for mean values or @c random for random samples.
 * @return A vector of weighted length of stay in each row in R_disease_prog. These
 * values are then summed by @c patient_id using @c data.table at the @c R level
 *  in the private member function @c IndivCtstm$sim_wlos. 
 ******************************************************************************/ 
// [[Rcpp::export]]
std::vector<double> C_indiv_ctstm_starting(Rcpp::DataFrame R_disease_prog,
                                           std::vector<int> strategy_idx,
                                           std::vector<int> patient_idx,
                                           Rcpp::Environment R_StateVal,
                                           double dr, std::string type
                                           ){
  hesim::ctstm::disease_prog disease_prog(R_disease_prog);
  bool time_reset = Rcpp::as<bool>(R_StateVal["time_reset"]);
  hesim::statmods::obs_index obs_index(hesim::statmods::get_id_object(R_StateVal));
  hesim::statevals stvals(R_StateVal);
  
  int N = disease_prog.sample_.size();
  std::vector<double> out(N);
  int time_index = 0;
  for (unsigned int i = 0; i < N; ++i){
    double time = disease_prog.time_start_[i];
    obs_index.set_strategy_index(strategy_idx[i]);
    obs_index.set_patient_index(patient_idx[i]);
    obs_index.set_health_index(disease_prog.from_[i]);
    obs_index.set_time_index(0);
    
    if (!time_reset){
      if (i > 0 && 
          (disease_prog.patient_id_[i] != disease_prog.patient_id_[i-1] ||
          disease_prog.sample_[i] != disease_prog.sample_[i-1])){
          time_index = 0;
      }
      while (obs_index.get_time_stop() < time){
        ++time_index;
        obs_index.set_time_index(time_index);
      }
    }
    double yhat = stvals.sim(disease_prog.sample_[i], obs_index(), type);
    out[i] = exp(-dr * time) * yhat;
  }
  return out;
}

/***************************************************************************//** 
 * @ingroup ctstm
 * Simulate length of stay given simulated disease progression 
 * (i.e., a path through a multi-state model) from an individual-level model. 
 * @param R_disease_prog An R object of simulating disease progression generated
 * using C_ctstm_sim_disease.
 * @param strategy_idx The strategy index starting at 0.
 * @param patient_idx The patient index starting at 0.
 * @param dr The discount rate.
 * @return A vector of weighted length of stay in each row in R_disease_prog. These
 * values are then summed by @c patient_id using @c data.table at the @c R level
 *  in the private member function @c IndivCtstm$sim_wlos. 
 ******************************************************************************/ 
// [[Rcpp::export]]
std::vector<double> C_indiv_ctstm_los(Rcpp::DataFrame R_disease_prog,
                                       std::vector<int> strategy_idx,
                                       std::vector<int> patient_idx,
                                       double dr){
  hesim::ctstm::disease_prog disease_prog(R_disease_prog);
  
  int N = disease_prog.sample_.size();
  std::vector<double> los(N);
  for (int i = 0; i < N; ++i){
    double time = disease_prog.time_stop_[i] -  disease_prog.time_start_[i];
    los[i] = hesim::pv(1, dr,
                        disease_prog.time_start_[i],
                        disease_prog.time_start_[i] + time);
  }
  return los;
}



