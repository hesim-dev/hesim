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
  Rcpp::List R_input_mats = Rcpp::as<Rcpp::List>(R_StateVals["input_mats"]);
  hesim::statmods::obs_index obs_index(R_input_mats);

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
            int obs = obs_index(k, 0, i, h, time_index);
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
