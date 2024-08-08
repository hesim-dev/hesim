#include <boost/numeric/odeint.hpp>
#include <RcppArmadillo.h>
#include <hesim.h>
#include <hesim/ctstm/ctstm.h>
#include <hesim/statevals.h>
#include <hesim/dtstm.h>
#include <hesim/utils.h> // member_of()

// Boilerplate code to ensure that arma plays nicely with boost::numeric::odeint
namespace boost {
  namespace numeric {
    namespace odeint {
      template <>
      struct is_resizeable<arma::vec>
      {
	typedef boost::true_type type;
	const static bool value = type::value;
      };
      template <>
      struct same_size_impl<arma::vec, arma::vec>
      {
	static bool same_size(const arma::vec& x, const arma::vec& y)
	{
	  return x.size() == y.size();
	}
      };
      template<>
      struct resize_impl<arma::vec, arma::vec>
      {
	static void resize(arma::vec& v1, const arma::vec& v2)
	{
	  v1.resize(v2.size());
	}
      };
      template <>
      struct is_resizeable<arma::mat>
      {
	typedef boost::true_type type;
	const static bool value = type::value;
      };
      template <>
      struct same_size_impl<arma::mat, arma::mat>
      {
	static bool same_size(const arma::mat& x, const arma::mat& y)
	{
	  return x.n_rows == y.n_rows && x.n_cols == y.n_cols;
	}
      };
      template<>
      struct resize_impl<arma::mat, arma::mat>
      {
	static void resize(arma::mat& m1, const arma::mat& m2)
	{
	  m1.resize(m2.n_rows, m2.n_cols);
	}
      };
    }
  }
} // namespace boost::numeric::odeint

/***************************************************************************//** 
 * @ingroup ctstm
 * Simulate disease progression (i.e., a path through a multi-state model),
 * together with costs and QALYs.
 * This function is exported to @c R and used in @c CohortCtstm$sim().
 * @param R_CtstmTrans An R object of class @c CohortCtstmTrans.
 * @param R_CostsStateVals An R List of @c StateVal objects for costs.
 * @param R_QALYsStateVal An R object of class @c StateVal for QALYs.
 * @param live_states An integer vector the for the live states.
 * @param start_state The starting health state for each patient,
 * (TODO: allow for probabilities in the initial states).
 * @param start_age The starting age of each patient in the simulation.
 * @param times A vector of times to be simulated.
 * @param clock A string to specify the clock; allows for "forward" and "mixt".
 * @param transition_types A vector of integers for the transitions when clock="mixt";
 * 0=Time on study, 1=Age.
 * @param progress An integer for how to report progress; defaults to zero, 
 * which is no reporting.
 * @param dr_qalys A double for the discount rate for QALYs (NB: this is not annualised).
 * @param dr_costs A double for the discount rate for costs (NB: this is not annualised).
 * @param type A string for how costs and QALYs are calculated; defaults to "predict".
 * @param eps A double for how to represent zero time; defaults to 1e-100.
 * @return An R List with elements "stateprobs", which is a data frame of the same format 
 * as stateprobs_out, and "ev", which is a data frame of the same format as ev_out.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::List C_cohort_ctstm_sim(Rcpp::Environment R_CtstmTrans,
			      Rcpp::List R_CostsStateVals,
			      Rcpp::Environment R_QALYsStateVal,
			      std::vector<int> live_states,
			      std::vector<int> start_state, // by patient
			      std::vector<double> start_age, // by patient
			      std::vector<double> times, // common
			      std::string clock,
			      std::vector<int> transition_types,
			      int progress = 0,
			      double dr_qalys = 0.0,
			      double dr_costs = 0.0,
			      std::string type="predict",
			      double zero_tol = 1e-100,
			      double abs_tol = 1e-6,
			      double rel_tol = 1e-6) {
  using namespace boost::numeric::odeint;
  typedef arma::vec state_type;
  enum TransitionType {tt_reset, tt_time, tt_age}; // tt_reset is *not* allowed here
  // Initialize
  std::unique_ptr<hesim::ctstm::transmod> transmod = hesim::ctstm::transmod::create(R_CtstmTrans);
  size_t n_costs = R_CostsStateVals.size();
  std::vector<hesim::statevals> costs_lookup;
  std::vector<hesim::statmods::obs_index> obs_index_costs;
  for (size_t cost_=0; cost_<n_costs; ++cost_) {
    costs_lookup.push_back(hesim::statevals(Rcpp::as<Rcpp::Environment>(R_CostsStateVals[cost_])));
    obs_index_costs.push_back(hesim::statmods::obs_index(hesim::statmods::get_id_object(Rcpp::as<Rcpp::Environment>(R_CostsStateVals[cost_]))));
  }
  // use unique_ptr because statevals and obs_index can be null and do not have null constructors
  std::unique_ptr<hesim::statevals> qalys_lookup;
  std::unique_ptr<hesim::statmods::obs_index> obs_index_qalys;
  bool valid_R_QALYsStateVal = R_QALYsStateVal.exists("method");
  if (valid_R_QALYsStateVal) {
    qalys_lookup = std::unique_ptr<hesim::statevals>(new hesim::statevals(R_QALYsStateVal));
    auto id_object = hesim::statmods::get_id_object(R_QALYsStateVal);
    obs_index_qalys =
      std::unique_ptr<hesim::statmods::obs_index>(new hesim::statmods::obs_index(id_object));
  }
  size_t n_samples = transmod->get_n_samples();
  size_t n_strategies = transmod->get_n_strategies(); 
  size_t n_patients = transmod->get_n_patients(); // actually, the number of patient profiles
  size_t n_times = times.size();
  size_t n_states = transmod->trans_mat_.n_states_; // NB: includes death states
  // Rprintf("n_states=%i\n",n_states);
  size_t N = n_samples * 
          n_strategies *  
          n_patients *
          n_states *
          n_times;
  hesim::stateprobs_out out(N);
  size_t N2 = n_samples * 
          n_strategies *  
          n_patients *
          n_states *
    (1+(valid_R_QALYsStateVal ? 1 : 0)+n_costs); // LYs, QALYs, cost categories
  size_t n_trans = transmod->trans_mat_.n_trans_;
  arma::uvec from(n_trans);
  arma::uvec to(n_trans);
  for (size_t state_id : live_states) {
    std::vector<int> trans_ids = transmod->trans_mat_.trans_id(state_id);
    std::vector<int> tos = transmod->trans_mat_.to(state_id);
    size_t n_trans_state = trans_ids.size();
    for (size_t trans_ = 0; trans_ < n_trans_state; ++trans_){ // NB: within a state
      size_t trans_id = trans_ids[trans_];
      from(trans_id) = state_id;
      to(trans_id) = tos[trans_];
    }
  }
  hesim::ev_out out2(N2);
  size_t counter = 0, counter2 = 0;
  // Rprintf("Main loop\n");
  // main loop
  for (size_t s = 0; s < n_samples; ++s){
    if (progress > 0){
      if ((s + 1) % progress == 0){ // R-based indexing
        Rcpp::Rcout << "sample = " << s + 1 << std::endl;  
      }
    }
    for (size_t k = 0; k < n_strategies; ++k){
      transmod->obs_index_.set_strategy_index(k);
        for (size_t i = 0; i < n_patients; ++i){
          transmod->obs_index_.set_patient_index(i);
	  // allow for point mass distributions using clock-forward
	  std::vector<bool> is_point_mass(n_trans);
	  std::vector<double> point_mass_times(n_trans);
	  std::vector<double> full_times(times);
	  // apologies for the ugly dynamic_cast's
	  hesim::ctstm::mstate_list* tmp_transmod =
	    dynamic_cast<hesim::ctstm::mstate_list*>(transmod.get());
	  for (size_t trans_ = 0; trans_ < n_trans; ++trans_) {
	    if (tmp_transmod) {
	      hesim::stats::point_mass* tmp_point_mass =
		dynamic_cast<hesim::stats::point_mass*>(tmp_transmod->get_dist(trans_));
	      if (tmp_point_mass) {
		is_point_mass[trans_] = true;
		point_mass_times[trans_] = tmp_point_mass->random(); // actually, it is fixed:)
		full_times.push_back(point_mass_times[trans_]);
	      } else {
		is_point_mass[trans_] = false;
		point_mass_times[trans_] = 0.0;
	      }
	    } else {
	      is_point_mass[trans_] = false;
	      point_mass_times[trans_] = 0.0;
	    }
	  }
	  // sort and unique for full_times
	  std::set<double> set_full_times(full_times.begin(), full_times.end());
	  full_times.assign(set_full_times.begin(), set_full_times.end());
	  size_t n_full_times = full_times.size();
	  arma::vec p0 = arma::zeros(n_states*(3+R_CostsStateVals.size()));
	  p0[start_state[i]] = 1.0;
	  arma::mat report(n_times,p0.size());
	  report.row(0) = p0.t();
	  auto stepper = make_dense_output(rel_tol, abs_tol, runge_kutta_dopri5<state_type>());
	  size_t j_counter=1;
	  for (size_t j=1; j<n_full_times; j++) {
	    // special case: point mass
	    for (size_t trans_ = 0; trans_ < n_trans; ++trans_) {
	      if (is_point_mass[trans_] && full_times[j-1]==point_mass_times[trans_]) {
		for (size_t cost_=0; cost_<n_costs; ++cost_) {
		  if (costs_lookup[cost_].method_ == "transition") {
		    size_t t_index = hesim::hesim_bound(full_times[j-1], obs_index_costs[cost_].time_start_);
		    size_t obs = obs_index_costs[cost_](k, // strategy
						     i, // patient
						     trans_, // assumes transition
						     t_index);
		    p0(from[trans_] + (3+cost_)*n_states) +=
		      p0(from[trans_])*costs_lookup[cost_].sim(s, obs, type);
		  }
		}
		// assumes transitions done in order if the times are the same
		p0(to[trans_]) += p0(from[trans_]);
		p0(from[trans_]) = 0.0;
	      }
	    }
	    size_t n = integrate_adaptive(stepper,
		       [&](const state_type &Y , state_type &dYdt, const double t)
		       {
			 arma::vec rates(n_trans);
			 arma::vec utilities(n_states);
			 arma::mat costs_by_state(n_states,n_costs);
			 arma::mat costs_by_transition(n_trans,n_costs);
			 double age = std::max(zero_tol, start_age[i]+t);
			 double tstar = std::max(zero_tol,t);
			 // Rprintf("Rate calculations\n");
			 // rate calculations
			 for (size_t trans_ = 0; trans_ < n_trans; ++trans_) {
			   // special case: point_mass
			   if (is_point_mass[trans_]) {
			       rates[trans_] = 0.0;
			   } else { // general case
			     if (clock == "forward"){
			       rates[trans_] = transmod->summary(trans_, s, {tstar}, "hazard")[0];
			     } else { // clock == "mixt"
			       if (transition_types[trans_] == tt_age){
				 rates[trans_] = transmod->summary(trans_, s, {age}, "hazard")[0];
			       } else { // transition_types[trans_] == tt_time (tt_reset is *not* currently available)
				 rates[trans_] = transmod->summary(trans_, s, {tstar}, "hazard")[0];
			       }
			     }
			   } // end general case
			 } // end loop over transitions
			 // Rprintf("Utility input calculations\n");
			 // utility input calculations
			 if (valid_R_QALYsStateVal) {
			   size_t t_index = hesim::hesim_bound(t, obs_index_qalys->time_start_);
			   for (size_t state_ : live_states) {
			     size_t obs = (*obs_index_qalys)(k, // strategy
							  i, // patient
							  state_,
							  t_index);
			     // Rprintf("Utility: k=%i, i=%i, state_=%i, t_index=%i, s=%i, obs=%i\n",
			     //     k, i, state_, t_index, s, obs);
			     utilities(state_) = qalys_lookup->sim(s, obs, type);
			     // Rprintf("Utility done\n");
			   } // end loop over states
			 }
			 // Rprintf("Cost input calculations\n");
			 // cost input calculations
			 for (size_t cost_=0; cost_<n_costs; ++cost_) {
			   size_t t_index = hesim::hesim_bound(t, obs_index_costs[cost_].time_start_);
			   if (costs_lookup[cost_].method_ == "wlos") {
			     for (size_t state_ : live_states) {
			       size_t obs = obs_index_costs[cost_](k, // strategy
								i, // patient
								state_,
								t_index);
			       costs_by_state(state_, cost_) += costs_lookup[cost_].sim(s, obs, type);
			     } // end loop over live states
			   } else if (costs_lookup[cost_].method_ == "starting") {
			     for (size_t trans_=0; trans_<n_trans; ++trans_) {
			       if (hesim::member_of(to[trans_], live_states)) {
				 size_t obs = obs_index_costs[cost_](k, // strategy
								  i, // patient
								  to(trans_),
								  t_index);
				 costs_by_transition(trans_, cost_) = costs_lookup[cost_].sim(s, obs, type);
			       }
			     } // end loop over transitions
			   } else { // costs_lookup[cost_].method_ == "transition"
			     for (size_t trans_=0; trans_<n_trans; ++trans_) {
			       size_t obs = obs_index_costs[cost_](k, // strategy
								i, // patient
								trans_, // assumes transition
								t_index);
			       costs_by_transition(trans_, cost_) = costs_lookup[cost_].sim(s, obs, type); // by trans_
			     } // end loop over transitions
			   } // end case: "transition"
			 } // end loop over costs
			 // Rprintf("Transition intensity calculations\n");
			 dYdt = dYdt*0.0;
			 arma::vec delta = Y(from) % rates;
			 dYdt(to) += delta;
			 dYdt(from) -= delta;
			 double drr_utilities = std::exp(-dr_qalys*t);
			 double drr_costs = std::exp(-dr_costs*t);
			 dYdt(arma::span(n_states,2*n_states-1)) = Y(arma::span(0,n_states-1));
			 // Rprintf("Discounted utilities\n");
			 // discounted utilities
			 dYdt(arma::span(2*n_states,3*n_states-1)) = utilities % Y(arma::span(0,n_states-1))*drr_utilities;
			 for (size_t cost_=0; cost_<n_costs; ++cost_) {
			   if (costs_lookup[cost_].method_ == "wlos")
			     dYdt(arma::span((3+cost_)*n_states,(4+cost_)*n_states-1)) += costs_by_state.col(cost_) % Y(arma::span(0,n_states-1))*drr_costs;
			   else if (costs_lookup[cost_].method_ == "starting") {
			     for (size_t trans_=0; trans_<n_trans; ++trans_)
			       dYdt(to[trans_] + (3+cost_)*n_states) += costs_by_transition(trans_,cost_) * Y(from[trans_]) * rates(trans_) *drr_costs;
			   } else { // costs_lookup[cost_].method_ == "transition"
			     for (size_t trans_=0; trans_<n_trans; ++trans_)
			       dYdt(from[trans_] + (3+cost_)*n_states) += // costs arbitrarily assigned to the state of origin (from state)
				 costs_by_transition(trans_,cost_) * Y(from(trans_)) * rates(trans_) * drr_costs;
			   }
			 }
		       },
					  p0,
					  full_times[j-1],
					  full_times[j],
					  full_times[j]-full_times[j-1]);
	    if (hesim::member_of(full_times[j], times)) {
	      report.row(j_counter) = p0.t();
	      j_counter++;
	    }
	  } // end j full_times_ loop
	  for (size_t h = 0; h < n_states; ++h){
	    for (size_t ti = 0; ti < n_times; ++ti){
	      out.sample_[counter] = s;
	      out.strategy_id_[counter] = k;
	      out.patient_id_[counter] = i;
	      out.grp_id_[counter] = transmod->obs_index_.get_grp_id();
	      out.patient_wt_[counter] = transmod->obs_index_.get_patient_wt();
	      out.state_id_[counter] = h;
	      out.t_[counter] = times[ti];
	      out.prob_[counter] = report(ti, h);
	      ++counter;                   
	    } // end cycles loop
	    // report for life-years
	    auto update_ev_out = [&](hesim::ev_out & object, size_t counter) {
	      object.sample_[counter] = s;
	      object.strategy_id_[counter] = k;
	      object.patient_id_[counter] = i;
	      object.grp_id_[counter] = transmod->obs_index_.get_grp_id();
	      object.patient_wt_[counter] = transmod->obs_index_.get_patient_wt();
	      object.state_id_[counter] = h;
	    };
	    update_ev_out(out2,counter2);
	    out2.dr_[counter2] = 0.0; // assume no discounting for life-years
	    out2.outcome_[counter2] = "ly";
	    out2.value_[counter2] = report(n_times-1,n_states+h);
	    ++counter2;
	    // report for qalys
	    if (valid_R_QALYsStateVal) {
	      update_ev_out(out2,counter2);
	      out2.dr_[counter2] = dr_qalys;
	      out2.outcome_[counter2] = "qaly";
	      out2.value_[counter2] = report(n_times-1,2*n_states+h);
	      ++counter2;
	    }
	    // report for costs
	    for (size_t cost_=0; cost_<n_costs; ++cost_) {
	      update_ev_out(out2,counter2);
	      out2.dr_[counter2] = dr_costs;
	      out2.outcome_[counter2] = "Category " + std::to_string(cost_+1);
	      out2.value_[counter2] = report(n_times-1,(3+cost_)*n_states+h);
	      ++counter2;
	    } // end cost category loop
	  } // end state loop
	} // end patient loop
    } // end strategy loop
  } // end parameter sampling loop
  // Return
  using namespace Rcpp;
  return(List::create(_["stateprobs"]=out.create_R_data_frame(),
		      _["ev"]=out2.create_R_data_frame()));
}

// test example that does not use hesim
class CohortCtstmTestODE {
public:
  typedef arma::vec state_type;
  size_t n_states;
  double discount_rate;
  CohortCtstmTestODE(double discount_rate = 0.0)
    : n_states(4), discount_rate(discount_rate) { }
  inline double hweibull(const double t, const double shape, const double scale) {
    return R::dweibull(std::max(1e-100,t),shape,scale,0)/
      R::pweibull(std::max(1e-100,t),shape,scale,0,0);
  }
  void operator() (const state_type &Y , state_type &dYdt, const double t)
  {
    run_step(Y, dYdt, t);
  }
  void run_step(const state_type &Y , state_type &dYdt, const double t)
  {
    using namespace arma;
    // conv_to<uvec>::from(an_ivec)
    // vec rates = {0.1, 0.2, 0.3};
    vec rates = {hweibull(t,0.8,1.0), hweibull(t,0.8,1.0),
      hweibull(t,0.8,2.0), hweibull(t,0.8,3.0)};
    vec utilities = {0.9, 0.8, 0.7, 0.0};
    vec costs_per_unit_time = {10000.0, 20000.0, 30000.0, 0.0}; // by state
    vec costs_per_transition = {10.0, 15.0, 20.0, 30.0}; // by transition
    vec starting_costs = {0.0, 2000.0, 3000.0, 0.0}; // by state
    uvec from = {0, 1, 1, 2};
    uvec to = {1, 0, 2, 3};
    dYdt = dYdt*0.0;
    // state occupation/transition probabilities
    vec delta = Y(from) % rates;
    dYdt(to) += delta;
    dYdt(from) -= delta;
    // double dr = 1.0/pow(1.0+discount_rate,t);
    double dr = std::exp(-discount_rate*t);
    // undiscounted life-years
    dYdt(span(n_states,2*n_states-1)) = Y(span(0,n_states-1));
    // discounted utilities
    dYdt(span(2*n_states,3*n_states-1)) = utilities % Y(span(0,n_states-1))*dr;
    // discounted costs
    dYdt(span(3*n_states,4*n_states-1)) = costs_per_unit_time % Y(span(0,n_states-1))*dr;
    // discounted transition costs
    for (size_t j=0; j<4; ++j)
      dYdt(from[j] + 4*n_states) += costs_per_transition[j] * Y(from[j]) * rates(j) * dr;
    // discounted starting costs
    for (size_t j=0; j<4; ++j)
      dYdt(to[j] + 5*n_states) += starting_costs[to[j]] * Y(from[j]) * rates(j) * dr;
  }
  Rcpp::List run(size_t start_state = 0, arma::vec times = {0.0, 1.0, 2.0, 10.0},
		 double discount_rate = 0.0) {
    this->discount_rate = discount_rate;
    using namespace boost::numeric::odeint;
    size_t n_times = times.size();
    arma::vec p0(n_states*6); // probs, lys, qalys, 3*costs
    p0(start_state) = 1.0;
    // combine the results
    arma::mat combined(n_times,p0.size());
    combined.row(0) = p0.t();
    auto stepper = make_dense_output(1.0e-6, 1.0e-6, runge_kutta_dopri5<state_type>());
    for (size_t step_=1; step_<n_times; step_++) {
      size_t n = integrate_adaptive(stepper,
				 [this](const state_type &Y , state_type &dYdt, const double t) {
				   this->run_step(Y,dYdt,t);
				 },
				 p0,
				 times[step_-1],
				 times[step_],
				 times[step_]-times[step_-1]);
      combined.row(step_) = p0.t();
    }
    size_t counter=0;
    hesim::stateprobs_out out(n_times * n_states);
    for (size_t h = 0; h < n_states; ++h){
      for (size_t ti = 0; ti < n_times; ++ti){
	out.sample_[counter] = 0; // C-based indexing
	out.strategy_id_[counter] = 0;
	out.patient_id_[counter] = 0;
	out.grp_id_[counter] = 1; // ???
	out.patient_wt_[counter] = 1.0;
	out.state_id_[counter] = h;
	out.t_[counter] = times[ti];
	out.prob_[counter] = combined(ti, h);
	++counter;                   
      } // end cycles loop
    } // end 
    hesim::ev_out out2(n_states * 5);
    size_t counter2 = 0;
    for (size_t h = 0; h < n_states; ++h) {
      auto update_ev_out = [&](hesim::ev_out & object,size_t counter) {
	object.sample_[counter] = 0;
	object.strategy_id_[counter] = 0;
	object.patient_id_[counter] = 0;
	object.grp_id_[counter] = 1;
	object.patient_wt_[counter] = 1.0;
	object.state_id_[counter] = h;
	object.dr_[counter] = discount_rate;
      };
      update_ev_out(out2,counter2);
      out2.dr_[counter2] = 0.0;
      out2.outcome_[counter2] = "ly";
      out2.value_[counter2] = combined(n_times-1,1*n_states+h);
      ++counter2;
      update_ev_out(out2,counter2);
      out2.outcome_[counter2] = "qaly";
      out2.value_[counter2] = combined(n_times-1,2*n_states+h);
      ++counter2;
      update_ev_out(out2,counter2);
      out2.outcome_[counter2] = "Category 1";
      out2.value_[counter2] = combined(n_times-1,3*n_states+h);
      ++counter2;
      update_ev_out(out2,counter2);
      out2.outcome_[counter2] = "Category 2";
      out2.value_[counter2] = combined(n_times-1,4*n_states+h);
      ++counter2;
      update_ev_out(out2,counter2);
      out2.outcome_[counter2] = "Category 3";
      out2.value_[counter2] = combined(n_times-1,5*n_states+h);
      ++counter2;
    }
    using namespace Rcpp;
    return List::create(_["stateprobs_"]=out.create_R_data_frame(),
			_["ev_"]=out2.create_R_data_frame());
  }
};
// [[Rcpp::export]]
Rcpp::List runCohortCtstmTestODE(size_t start_state, arma::vec times,
				 double discount_rate = 0.0) {
  return CohortCtstmTestODE().run(start_state,times,discount_rate);
}

/***************************************************************************//** 
 * @ingroup ctstm
 * Simulate disease progression (i.e., a path through a semi-Markov multi-state model),
 * together with costs and QALYs.
 * This function is exported to @c R and used in @c CohortCtstm$sim().
 * @param R_CtstmTrans An R object of class @c CohortCtstmTrans.
 * @param R_CostsStateVals An R List of @c StateVal objects for costs.
 * @param R_QALYsStateVal An R object of class @c StateVal for QALYs.
 * @param live_states An integer vector the for the live states.
 * @param start_state The starting health state for each patient,
 * (TODO: allow for probabilities in the initial states).
 * @param start_age The starting age of each patient in the simulation.
 * @param max_time Double for the maximum time to be simulated.
 * @param clock A string to specify the clock; allows for "forward" and "mixt".
 * @param transition_types A vector of integers for the transitions when clock="mixt";
 * 0=Time on study, 1=Age.
 * @param delta_time Double for the time interval used to discretise time for the semi-Markov
 * models.
 * @param progress An integer for how to report progress; defaults to zero, 
 * which is no reporting.
 * @param dr_qalys A double for the discount rate for QALYs (NB: this is not annualised).
 * @param dr_costs A double for the discount rate for costs (NB: this is not annualised).
 * @param type A string for how costs and QALYs are calculated; defaults to "predict".
 * @param eps A double for how to represent zero time; defaults to 1e-100.
 * @param debug A boolean whether to print out debugging information
 * @return An R List with elements "stateprobs", which is a data frame of the same format 
 * as stateprobs_out, and "ev", which is a data frame of the same format as ev_out.
 ******************************************************************************/ 
// [[Rcpp::export]]
Rcpp::List C_cohort_ctstm_sim2(Rcpp::Environment R_CtstmTrans,
			       Rcpp::List R_CostsStateVals,
			       Rcpp::Environment R_QALYsStateVal,
			       std::vector<int> live_states,
			       std::vector<int> start_state, // by patient
			       std::vector<double> start_age, // by patient
			       double max_time,
			       std::string clock,
			       std::vector<int> transition_types,
			       double delta_time = 0.1, // new
			       int progress = 0,
			       double dr_qalys = 0.0,
			       double dr_costs = 0.0,
			       std::string type="predict",
			       double zero_tol = 1e-100,
			       double abs_tol = 1e-6,
			       double rel_tol = 1e-6,
			       bool debug = false) {
  using namespace boost::numeric::odeint;
  typedef arma::vec state_type;
  enum TransitionType {tt_reset, tt_time, tt_age};
  // Initialize
  std::unique_ptr<hesim::ctstm::transmod> transmod = hesim::ctstm::transmod::create(R_CtstmTrans);
  size_t n_costs = R_CostsStateVals.size();
  std::vector<hesim::statevals> costs_lookup;
  std::vector<hesim::statmods::obs_index> obs_index_costs;
  std::vector<bool> costs_time_reset;
  for (size_t cost_=0; cost_<n_costs; ++cost_) {
    costs_time_reset.push_back(Rcpp::as<bool>(Rcpp::as<Rcpp::Environment>(R_CostsStateVals[cost_])["time_reset"]));
    costs_lookup.push_back(hesim::statevals(Rcpp::as<Rcpp::Environment>(R_CostsStateVals[cost_])));
    obs_index_costs.push_back(hesim::statmods::obs_index(hesim::statmods::get_id_object(Rcpp::as<Rcpp::Environment>(R_CostsStateVals[cost_]))));
  }
  // use unique_ptr because statevals and obs_index can be null and do not have null constructors
  std::unique_ptr<hesim::statevals> qalys_lookup;
  std::unique_ptr<hesim::statmods::obs_index> obs_index_qalys;
  bool valid_R_QALYsStateVal = R_QALYsStateVal.exists("method");
  bool qalys_time_reset = false;
  if (valid_R_QALYsStateVal) {
    qalys_time_reset = Rcpp::as<bool>(R_QALYsStateVal["time_reset"]);
    qalys_lookup = std::unique_ptr<hesim::statevals>(new hesim::statevals(R_QALYsStateVal));
    auto id_object = hesim::statmods::get_id_object(R_QALYsStateVal);
    obs_index_qalys =
      std::unique_ptr<hesim::statmods::obs_index>(new hesim::statmods::obs_index(id_object));
  }
  size_t n_samples = transmod->get_n_samples();
  size_t n_strategies = transmod->get_n_strategies(); 
  size_t n_patients = transmod->get_n_patients(); // actually, the number of patient profiles
  size_t n_times = ceil(max_time / delta_time);
  size_t n_states = transmod->trans_mat_.n_states_; // NB: includes death states
  size_t n_reps = 2+(valid_R_QALYsStateVal ? 1 : 0)+n_costs;
  if (debug) Rprintf("n_states=%zu\n",n_states);
  size_t N = n_samples * 
    n_strategies *  
    n_patients *
    n_states *
    n_times;
  hesim::stateprobs_out out(N);
  size_t N2 = n_samples * 
    n_strategies *  
    n_patients *
    n_states *
    (1+(valid_R_QALYsStateVal ? 1 : 0)+n_costs); // LYs, QALYs, cost categories
  size_t n_trans = transmod->trans_mat_.n_trans_;
  arma::uvec from(n_trans);
  arma::uvec to(n_trans);
  for (size_t state_id : live_states) {
    std::vector<int> trans_ids = transmod->trans_mat_.trans_id(state_id);
    std::vector<int> tos = transmod->trans_mat_.to(state_id);
    size_t n_trans_state = trans_ids.size();
    for (size_t trans_ = 0; trans_ < n_trans_state; ++trans_){ // NB: within a state
      size_t trans_id = trans_ids[trans_];
      from(trans_id) = state_id;
      to(trans_id) = tos[trans_];
    }
  }
  hesim::ev_out out2(N2);
  size_t counter = 0, counter2 = 0;
  if (debug) Rprintf("Main loop\n");
  // main loop
  for (size_t s = 0; s < n_samples; ++s){
    if (progress > 0){
      if ((s + 1) % progress == 0){ // R-based indexing
        Rcpp::Rcout << "sample = " << s + 1 << std::endl;  
      }
    }
    for (size_t k = 0; k < n_strategies; ++k){
      transmod->obs_index_.set_strategy_index(k);
      for (size_t i = 0; i < n_patients; ++i){
	transmod->obs_index_.set_patient_index(i);
	// NOTE: can utilities and costs depend on duration -- or can we collapse for duration?
	size_t n_rows = n_states*(3+R_CostsStateVals.size());
	arma::vec p0 = arma::zeros(n_rows);
	p0(start_state[i]) = 1.0; // assumes starts with zero duration (TODO: generalise)
	arma::mat report(n_rows,n_times); // transposed cf. C_cohort_ctstm_sim()
	report.col(0) = p0;
	auto stepper = make_dense_output(rel_tol, abs_tol, runge_kutta_dopri5<state_type>());
	for (size_t j=1; j<n_times; j++) {
	  if (debug) Rprintf("j=%zu\n",j);
	  p0 = join_vert(p0, arma::zeros(n_rows)); // add zeros at the end for flows out
	  size_t n = integrate_adaptive(stepper,
					[&](const state_type &Y , state_type &dYdt, const double u)
					{
					  arma::mat y(Y.memptr(), Y.n_elem/(j+1), j); // NB: excludes the last column
					  arma::vec y0 = arma::zeros(n_rows); 
					  arma::mat dydt = y*0.0;
					  arma::vec dy0dt = y0*0.0;
					  arma::vec rates0(n_trans);
					  arma::mat rates(n_trans,j);
					  arma::mat utilities(n_states,j);
					  arma::cube costs_by_state(n_states,j,n_costs);
					  arma::cube costs_by_transition(n_trans,j,n_costs);
					  double t = std::max(zero_tol, (j-1)*delta_time+u);
					  double age = start_age[i]+t;
					  arma::vec duration = arma::linspace(0.0, (j-1)*delta_time, j)+u;
					  duration(0) = std::max(zero_tol, duration(0));
					  std::vector<double> vduration(duration.begin(), duration.end());
					  if (debug) Rprintf("Rate calculations\n");
					  for (size_t trans_ = 0; trans_ < n_trans; ++trans_) {
					    if (debug) Rprintf("trans_=%zu\n",trans_);
					    if (clock == "forward"){ // NB: this should never be reached
					      if (debug) Rprintf("clock forward\n");
					      rates.row(trans_) = rates.row(trans_)*0.0+transmod->summary(trans_, s, {t}, "hazard")[0];
					      rates0(trans_) = transmod->summary(trans_, s, {t-u/2}, "hazard")[0];
					    } else if (clock == "mixt") {
					      if (debug) Rprintf("clock mixt\n");
					      if (transition_types[trans_] == tt_age) {
						if (debug) Rprintf("tt_age\n");
						rates.row(trans_) = rates.row(trans_)*0.0+transmod->summary(trans_, s, {age}, "hazard")[0];
						rates0(trans_) = transmod->summary(trans_, s, {age-u/2}, "hazard")[0];
					      } else if (transition_types[trans_] == tt_time) {
						if (debug) Rprintf("tt_time\n");
						rates.row(trans_) = rates.row(trans_)*0.0+transmod->summary(trans_, s, {t}, "hazard")[0];
						rates0(trans_) = transmod->summary(trans_, s, {t-u/2}, "hazard")[0];
					      } else { // tt_reset
						if (debug) Rprintf("tt_reset\n");
						std::vector<double> vrow =
						  transmod->summary(trans_, s, vduration, "hazard");
						if (debug) Rprintf("vrow.size()=%zu, rates.n_rows=%i, rates.n_cols=%i\n",
								   vrow.size(), rates.n_rows, rates.n_cols);
						rates.row(trans_) = arma::rowvec(vrow.data(), vrow.size(), false);
						rates0(trans_) = transmod->summary(trans_, s, {u/2.0}, "hazard")[0];
					      }
					    } else { // clock == "reset"
					      if (debug) Rprintf("clock reset\n");
					      std::vector<double> vrow =
						transmod->summary(trans_, s, vduration, "hazard");
					      if (debug) Rprintf("vrow.size()=%zu, rates.n_rows=%i, rates.n_cols=%i\n",
						      vrow.size(), rates.n_rows, rates.n_cols);
					      rates.row(trans_) = arma::rowvec(vrow.data(), vrow.size(), false);
					      rates0(trans_) = transmod->summary(trans_, s, {u/2.0}, "hazard")[0];
					    }
					    if (debug) Rprintf("end loop trans_\n");
					  } // end loop over transitions
					  // utility input calculations
					  if (valid_R_QALYsStateVal) {
					    if (debug) Rprintf("Utility input calculations\n");
					    for (size_t d_index=0; d_index < duration.size(); ++d_index) {
					      size_t t_index = qalys_time_reset ?
						hesim::hesim_bound(duration[d_index]+u,obs_index_qalys->time_start_) :
						hesim::hesim_bound(t,obs_index_qalys->time_start_);
					      for (size_t state_ : live_states) {
						size_t obs = (*obs_index_qalys)(k, // strategy
									     i, // patient
									     state_,
									     t_index);
						// if (debug) Rprintf("Utility: k=%zu, i=%zu, state_=%zu, t_index=%zu, s=%zu, obs=%zu\n",
						// 		     k, i, state_, t_index, s, obs);
						utilities(state_,d_index) = qalys_lookup->sim(s, obs, type);
					      } // loop over state_
					    } // loop over d_index
					    if (debug) Rprintf("Utility done\n");
					  }
					  if (debug) Rprintf("Cost input calculations\n");
					  for (size_t cost_=0; cost_<n_costs; ++cost_) {
					    if (debug) Rprintf("Cost %zu\n", cost_);
					    size_t t_index = 0;
					    if (!costs_time_reset[cost_]) {
					      t_index = hesim::hesim_bound(t, obs_index_costs[cost_].time_start_);
					      if (costs_lookup[cost_].method_ == "wlos") {
						for (size_t state_ : live_states) {
						  size_t obs = obs_index_costs[cost_](k, // strategy
										   i, // patient
										   state_,
										   t_index);
						  costs_by_state(arma::span(state_,state_), arma::span::all, arma::span(cost_,cost_)) +=
						    costs_by_state(arma::span(state_,state_), arma::span::all, arma::span(cost_,cost_))*0.0 + costs_lookup[cost_].sim(s, obs, type);
						} // end loop over live states
					      }
					      else if (costs_lookup[cost_].method_ == "starting") {
						for (size_t trans_=0; trans_<n_trans; ++trans_) {
						  if (hesim::member_of(to[trans_], live_states)) {
						    size_t obs = obs_index_costs[cost_](k, // strategy
										     i, // patient
										     to(trans_),
										     t_index);
						    costs_by_transition(arma::span(trans_,trans_), arma::span::all, arma::span(cost_,cost_)) =
						      costs_by_transition(arma::span(trans_,trans_), arma::span::all, arma::span(cost_,cost_))*0.0 + costs_lookup[cost_].sim(s, obs, type);
						  } // end condition for live states
						} // loop over trans_
					      } // end case "starting"
					      else { // costs_lookup[cost_].method_ == "transition"
						for (size_t trans_=0; trans_<n_trans; ++trans_) {
						  size_t obs = obs_index_costs[cost_](k, // strategy
										   i, // patient
										   trans_, // assumes transition
										   t_index);
						  costs_by_transition(arma::span(trans_,trans_), arma::span::all, arma::span(cost_,cost_)) =
						    costs_by_transition(arma::span(trans_,trans_), arma::span::all, arma::span(cost_,cost_))*0.0 + costs_lookup[cost_].sim(s, obs, type); // by trans_
						} // end loop over trans_
					      }
					    } else { 
					      for (size_t d_index = 0; d_index < duration.size(); ++d_index) {
						t_index = hesim::hesim_bound(duration[d_index]+u,
									     obs_index_costs[cost_].time_start_);
						if (costs_lookup[cost_].method_ == "wlos") {
						  for (size_t state_ : live_states) {
						    size_t obs = obs_index_costs[cost_](k, // strategy
										     i, // patient
										     state_,
										     t_index);
						    costs_by_state(state_, d_index, cost_) += costs_lookup[cost_].sim(s, obs, type);
						  } // end loop over live states
						}
						else if (costs_lookup[cost_].method_ == "starting") {
						  for (size_t trans_=0; trans_<n_trans; ++trans_) {
						    if (hesim::member_of(to[trans_], live_states)) {
						      size_t obs = obs_index_costs[cost_](k, // strategy
										       i, // patient
										       to(trans_),
										       t_index);
						      costs_by_transition(trans_, d_index, cost_) =
							costs_lookup[cost_].sim(s, obs, type);
						    } // end condition for live states
						  } // loop over trans_
						} // end case "starting"
						else { // costs_lookup[cost_].method_ == "transition"
						  for (size_t trans_=0; trans_<n_trans; ++trans_) {
						    size_t obs = obs_index_costs[cost_](k, // strategy
										     i, // patient
										     trans_, // assumes transition
										     t_index);
						    costs_by_transition(trans_, d_index, cost_) =
						      costs_lookup[cost_].sim(s, obs, type); // by trans_
						  } // end loop over trans_
						}
					      } // end loop over d_index
					    }
					  } // end loop over costs
					  if (debug) Rprintf("Transition intensity calculations\n");
					  // flow from y
					  arma::uvec all_columns = arma::regspace<arma::uvec>(0, y.n_cols-1);
					  if (debug) Rprintf("flow_y\n");
					  arma::mat flow_y = y(from,all_columns) % rates;
					  if (debug) Rprintf("dydt\n");
					  dydt(from,all_columns) -= flow_y;
					  if (debug) Rprintf("dy0dt\n");
					  dy0dt(to) += sum(flow_y,1);
					  // flow from y0 to y0
					  if (debug) Rprintf("Transition intensity (y0)\n");
					  arma::vec flow_y0 = y0(from) % rates0;
					  dy0dt(from) -= flow_y0;
					  dy0dt(to) += flow_y0;
					  // discounting
					  double drr_utilities = std::exp(-dr_qalys*t);
					  double drr_costs = std::exp(-dr_costs*t);
					  if (debug) Rprintf("Length of stay calculations\n");
					  dydt(arma::span(n_states,2*n_states-1), arma::span::all) =
					    y(arma::span(0,n_states-1), arma::span::all);
					  dy0dt(arma::span(n_states,2*n_states-1)) = y0(arma::span(0,n_states-1));
					  if (debug) Rprintf("Discounted utilities\n");
					  dydt(arma::span(2*n_states,3*n_states-1),arma::span::all) =
					    y(arma::span(0,n_states-1), arma::span::all) % utilities * drr_utilities;
					  if (debug) Rprintf("Discounted costs (y0)\n");
					  dy0dt(arma::span(2*n_states,3*n_states-1)) =
					    y0(arma::span(0,n_states-1)) % utilities.col(0) * drr_utilities;
					  if (debug) Rprintf("Discounted costs\n");
					  for (size_t cost_=0; cost_<n_costs; ++cost_) {
					    if (debug) Rprintf("Costs %zu\n", cost_);
					    if (costs_lookup[cost_].method_ == "wlos") {
					      if (debug) Rprintf("wlos\n");
					      arma::mat costs = costs_by_state.slice(cost_);
					      dydt(arma::span((3+cost_)*n_states,(4+cost_)*n_states-1),arma::span::all) +=
						y(arma::span(0,n_states-1),arma::span::all) %
						costs * drr_costs;
					      if (debug) Rprintf("y0 calculation\n");
					      dy0dt(arma::span((3+cost_)*n_states,(4+cost_)*n_states-1)) +=
						y0(arma::span(0,n_states-1)) % costs.col(0) * drr_costs;
					    } else if (costs_lookup[cost_].method_ == "starting") {
					      if (debug) Rprintf("starting\n");
					      arma::mat costs = costs_by_transition.slice(cost_);
					      for (size_t trans_=0; trans_<n_trans; ++trans_) {
						dydt.row(to(trans_) + (3+cost_)*n_states) +=
						  costs.row(trans_) % y.row(from(trans_)) % rates.row(trans_) * drr_costs;
						if (debug) Rprintf("y0 calculation\n");
						dy0dt(to(trans_) + (3+cost_)*n_states) +=
						  costs(trans_,0) * y0(from(trans_)) * rates0(trans_) * drr_costs;
					      }
					    } else { // costs_lookup[cost_].method_ == "transition"
					      if (debug) Rprintf("transition\n");
					      arma::mat costs = costs_by_transition.slice(cost_);
					      for (size_t trans_=0; trans_<n_trans; ++trans_) {
						dydt(from(trans_) + (3+cost_)*n_states, arma::span::all) +=
						  // costs arbitrarily assigned to the state of origin (from state)
						  costs.row(trans_) %
						  y.row(from(trans_)) %
						  rates.row(trans_) * drr_costs;
						if (debug) Rprintf("y0 calculation\n");
						dy0dt(from(trans_) + (3+cost_)*n_states) +=
						  costs(trans_,0) * y0(from(trans_)) * rates0(trans_) * drr_costs;
					      }
					    }
					  }
					  if (debug) Rprintf("End of ODE processing\n");
					  if (debug) Rprintf("dydt.n_rows=%i, dydt.n_cols=%i, dy0dt.n_elem=%i\n",
						  dydt.n_rows, dydt.n_cols, dy0dt.n_elem);
					  arma::mat dy = join_horiz(dydt, dy0dt);
					  if (debug) Rprintf("create vector...\n");
					  arma::vec dYdt_copy(dy.memptr(), dy.n_cols * dy.n_rows);
					  if (debug) Rprintf("assign copy...\n");
					  dYdt = dYdt_copy;
					  if (debug) Rprintf("End of ODE\n");
					},
					p0,
					0.0,
					delta_time,
					delta_time);
	  if (debug) Rprintf("Calculate new p0\n");
	  if (debug) Rprintf("p0.n_elem=%i, n_rows=%zu\n", p0.n_elem, n_rows);
	  arma::vec w0a = 0.5 * arma::ones(n_states);
	  for (size_t w0_index=0; w0_index < n_states; ++w0_index)
	    if (p0(p0.n_elem-n_rows+w0_index) != 0.0)
	      w0a(w0_index) = 1.0 - p0(p0.n_elem-n_rows+n_states+w0_index)/p0(p0.n_elem-n_rows+w0_index)/delta_time;
	  if (debug) w0a.print("w0a=");
	  arma::vec w0 = w0a;
	  for (size_t w0_index=1; w0_index<n_reps; ++w0_index)
	    w0 = join_vert(w0, w0a);
	  if (j>1)
	    p0 = join_cols(join_cols(w0 % p0(arma::span(p0.n_elem-n_rows,p0.n_elem-1)),  // half of y0
				     (1-w0) %  p0(arma::span(p0.n_elem-n_rows,p0.n_elem-1)) + // half of y0 and the first column
				     p0(arma::span(0,n_rows-1))),
			   p0(arma::span(n_rows,p0.n_elem-n_rows-1))); // the remainder of p0
	  else // j==1
	    p0 = join_cols(w0 % p0(arma::span(n_rows,p0.n_elem-1)),  // half of y0
			   (1-w0) % p0(arma::span(n_rows,p0.n_elem-1)) + // half of y0 and the first column
			   p0(arma::span(0,n_rows-1)));
	  if (debug) Rprintf("Add to report\n");
	  report.col(j) = sum(arma::mat(p0.memptr(), n_rows, j+1),1);
	  if (debug) Rprintf("End of j loop\n");
	} // end j loop
	if (debug) Rprintf("Start of reporting\n");
	for (size_t h = 0; h < n_states; ++h){
	  if (debug) Rprintf("Report stateprobs\n");
	  for (size_t ti = 0; ti < n_times; ++ti){
	    out.sample_[counter] = s;
	    out.strategy_id_[counter] = k;
	    out.patient_id_[counter] = i;
	    out.grp_id_[counter] = transmod->obs_index_.get_grp_id();
	    out.patient_wt_[counter] = transmod->obs_index_.get_patient_wt();
	    out.state_id_[counter] = h;
	    out.t_[counter] = delta_time*ti;
	    out.prob_[counter] = report(h,ti);
	    ++counter;
	  } // end cycles loop
	    // report for life-years
	  auto update_ev_out = [&](hesim::ev_out & object,size_t counter) {
	    object.sample_[counter] = s;
	    object.strategy_id_[counter] = k;
	    object.patient_id_[counter] = i;
	    object.grp_id_[counter] = transmod->obs_index_.get_grp_id();
	    object.patient_wt_[counter] = transmod->obs_index_.get_patient_wt();
	    object.state_id_[counter] = h;
	  };
	  if (debug) Rprintf("Report for life-years\n");
	  update_ev_out(out2,counter2);
	  out2.dr_[counter2] = 0.0; // assume no discounting for life-years
	  out2.outcome_[counter2] = "ly";
	  out2.value_[counter2] = report(n_states+h,n_times-1);
	  ++counter2;
	  if (debug) Rprintf("Report for QALYs\n");
	  // report for qalys
	  if (valid_R_QALYsStateVal) {
	    update_ev_out(out2,counter2);
	    out2.dr_[counter2] = dr_qalys;
	    out2.outcome_[counter2] = "qaly";
	    out2.value_[counter2] = report(2*n_states+h,n_times-1);
	    ++counter2;
	  }
	  if (debug) Rprintf("Report for costs\n");
	  for (size_t cost_=0; cost_<n_costs; ++cost_) {
	    update_ev_out(out2,counter2);
	    out2.dr_[counter2] = dr_costs;
	    out2.outcome_[counter2] = "Category " + std::to_string(cost_+1);
	    out2.value_[counter2] = report((3+cost_)*n_states+h,n_times-1);
	    ++counter2;
	  } // end cost category loop
	} // end state loop
      } // end patient loop
    } // end strategy loop
  } // end parameter sampling loop
  // Return
  using namespace Rcpp;
  return(List::create(_["stateprobs"]=out.create_R_data_frame(),
		      _["ev"]=out2.create_R_data_frame()));
}
