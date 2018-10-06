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
  std::vector<int> strategy_id = Rcpp::as<std::vector<int> >(R_input_mats["strategy_id"]);
  std::vector<int> patient_id = Rcpp::as<std::vector<int> >(R_input_mats["patient_id"]);
  std::vector<int> state_id = Rcpp::as<std::vector<int> >(R_input_mats["state_id"]);
  
  int n_samples = statevals.statmod_->get_n_samples();
  int n_obs = state_id.size();
  
  // Storage
  int N = times.size() * n_samples * n_obs;
  std::vector<int> state_id_vec(N);
  std::vector<int> sample_vec(N);
  std::vector<int> strategy_id_vec(N);
  std::vector<int> patient_id_vec(N);
  std::vector<double> times_vec(N);
  std::vector<double> value_vec(N);
  
  // Main loop
  int index = 0;
  for (int t = 0 ; t < times.size(); ++t){
    for (int s = 0; s < n_samples; ++s){
      for (int i = 0; i < n_obs; ++i){
        state_id_vec[index] = state_id[i];
        sample_vec[index] = s;
        strategy_id_vec[index] = strategy_id[i];
        patient_id_vec[index] = patient_id[i];
        times_vec[index] = times[t];
        value_vec[index] = statevals.sim(s, i, type);
        ++index;
      } // end observation loop
    } // end random samples loop
  } // end time loop
  
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
