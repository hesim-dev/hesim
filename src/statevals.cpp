#include <hesim/statevals.h>

/***************************************************************************//** 
 * @ingroup statevals
 * Simulate state values over time.
 * An Rcpp exported function that simulates state values over time using
 * @c statevals::sim.
 * @param R_StateVals
 * @param times
 * @param type
 * @return N
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::DataFrame C_statevals_sim(Rcpp::Environment R_StateVals, 
                           std::vector<double> times,
                           std::string type){
  // Initialize
  hesim::statevals statevals(R_StateVals);
  hesim::statmods::obs_index obs_index(hesim::statmods::get_id_object(R_StateVals));

  // Storage
  int n_samples = statevals.statmod_->get_n_samples();
  int N = times.size() * n_samples * obs_index.n_strategies_ *
          obs_index.n_patients_ * obs_index.n_healthvals_;
  std::vector<int> state_id_vec(N);
  std::vector<int> sample_vec(N);
  std::vector<int> strategy_id_vec(N);
  std::vector<int> patient_id_vec(N);
  std::vector<double> times_vec(N);
  std::vector<double> value_vec(N);
  
  // Main loop
  int index = 0;
  for (int s = 0; s < n_samples; ++s){
    for (int k = 0; k < obs_index.n_strategies_;++k){
      for (int i = 0; i < obs_index.n_patients_;++i){
        for (int h = 0; h < obs_index.n_healthvals_; ++h){
          int time_index = 0;
          for (int t = 0; t < times.size(); ++t){
            int obs = obs_index(k, i, h, time_index);
            strategy_id_vec[index] = obs_index.get_strategy_id();
            sample_vec[index] = s;
            state_id_vec[index] = obs_index.get_health_id();
            patient_id_vec[index] = obs_index.get_patient_id();
            times_vec[index] = times[t];
            value_vec[index] = statevals.sim(s, obs, type);
            ++index;
            if (times[t] >= obs_index.get_time_stop()){
              ++time_index;
            }
          } // end time loop
        } // end health states loop
      } // end patients loop
    } // end strategy loop
  } // end random samples loop
  
// Return
  return Rcpp::DataFrame::create(
    Rcpp::_["sample"] = sample_vec,
    Rcpp::_["strategy_id"] = strategy_id_vec,
    Rcpp::_["patient_id"] = patient_id_vec,
    Rcpp::_["state_id"] = state_id_vec,    
    Rcpp::_["time"] = times_vec,
    Rcpp::_["value"] = value_vec,
    Rcpp::_["stringsAsFactors"] = false
  );  
}

/***************************************************************************//** 
 * @ingroup statevals
 * Simulate weighted length of stay from simulated health state probabilities.
 * This function is exported to @c R and used to simulate costs and utilities.
 * @param R_stateprobs Simulated state probabilities from @c R.
 * @param statevals A list of @c R objects of class @c StateVals.
 * @param dr Discount rate.
 * @param categories Categories with a given @p type. For QALYs, there is only one
 * category ("qalys"), but for costs this could consist of different cost categories
 * such as drug acquisition and administration costs, resource use costs, etc. 
 * @return An @c R data frame with columns equivalent to the data members in
 * hesim::wlos_out_out.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::DataFrame C_sim_wlos(Rcpp::DataFrame R_stateprobs,
                               Rcpp::List R_statevals, 
                               std::vector<double> dr,
                               std::vector<std::string> categories,
                               std::vector<double> times,
                               std::string method = "trapz"){
  hesim::wlos wlos(R_statevals);
  hesim::stateprobs_out stprobs(R_stateprobs);
  hesim::wlos_out out = wlos(stprobs, times, dr, categories, method);
  return out.create_R_data_frame();
}

