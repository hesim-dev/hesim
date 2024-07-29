#include <boost/numeric/odeint.hpp>
#include <RcppArmadillo.h>
#include <hesim.h>
#include <hesim/ctstm/ctstm.h>
#include <hesim/statevals.h>
#include <hesim/dtstm.h>

#include <algorithm> // count

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
    }
  }
} // namespace boost::numeric::odeint

#define CTSTM_OUT(OBJECT,COUNTER)	     \
  OBJECT.sample_[COUNTER] = s;		     \
  OBJECT.strategy_id_[COUNTER] = k;	     \
  OBJECT.patient_id_[COUNTER] = i;					\
  OBJECT.grp_id_[COUNTER] = transmod->obs_index_.get_grp_id();		\
  OBJECT.patient_wt_[COUNTER] = transmod->obs_index_.get_patient_wt();	\
  OBJECT.state_id_[COUNTER] = h;

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
  enum TransitionType {tt_time, tt_age};
  // Initialize
  std::unique_ptr<hesim::ctstm::transmod> transmod = hesim::ctstm::transmod::create(R_CtstmTrans);
  int n_costs = R_CostsStateVals.size();
  std::vector<hesim::statevals> costs_lookup;
  std::vector<hesim::statmods::obs_index> obs_index_costs;
  for (int cost_=0; cost_<n_costs; ++cost_) {
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
  int n_samples = transmod->get_n_samples();
  int n_strategies = transmod->get_n_strategies(); 
  int n_patients = transmod->get_n_patients(); // actually, the number of patient profiles
  int n_times = times.size();
  int n_states = transmod->trans_mat_.n_states_; // NB: includes death states
  // Rprintf("n_states=%i\n",n_states);
  int N = n_samples * 
          n_strategies *  
          n_patients *
          n_states *
          n_times;
  hesim::stateprobs_out out(N);
  int N2 = n_samples * 
          n_strategies *  
          n_patients *
          n_states *
    (1+(valid_R_QALYsStateVal ? 1 : 0)+n_costs); // LYs, QALYs, cost categories
  int n_trans = transmod->trans_mat_.n_trans_;
  arma::uvec from(n_trans);
  arma::uvec to(n_trans);
  for (int state_id : live_states) {
    std::vector<int> trans_ids = transmod->trans_mat_.trans_id(state_id);
    std::vector<int> tos = transmod->trans_mat_.to(state_id);
    int n_trans_state = trans_ids.size();
    for (int trans_ = 0; trans_ < n_trans_state; ++trans_){ // NB: within a state
      int trans_id = trans_ids[trans_];
      from(trans_id) = state_id;
      to(trans_id) = tos[trans_];
    }
  }
  hesim::ev_out out2(N2);
  int counter = 0, counter2 = 0;
  // Rprintf("Main loop\n");
  // main loop
  for (int s = 0; s < n_samples; ++s){
    if (progress > 0){
      if ((s + 1) % progress == 0){ // R-based indexing
        Rcpp::Rcout << "sample = " << s + 1 << std::endl;  
      }
    }
    for (int k = 0; k < n_strategies; ++k){
      transmod->obs_index_.set_strategy_index(k);
        for (int i = 0; i < n_patients; ++i){
          transmod->obs_index_.set_patient_index(i);
	  arma::vec p0 = arma::zeros(n_states*(3+R_CostsStateVals.size()));
	  p0[start_state[i]] = 1.0;
	  arma::mat report(n_times,p0.size());
	  report.row(0) = p0.t();
	  auto stepper = make_dense_output(rel_tol, abs_tol, runge_kutta_dopri5<state_type>());
	  for (size_t j=1; j<n_times; j++) {
	    size_t n = integrate_adaptive(stepper,
		       [&](const state_type &Y , state_type &dYdt, const double t)
		       {
			 arma::vec rates(n_trans);
			 arma::vec utilities(n_states);
			 arma::mat costs_by_state(n_states,n_costs);
			 arma::mat costs_by_transition(n_trans,n_costs);
			 double age = start_age[i]+t;
			 double tstar = std::max(zero_tol,t);
			 // Rprintf("Rate calculations\n");
			 // rate calculations
			 for (int state_ = 0; state_<n_states; ++state_) {
			   for (int trans_ = 0; trans_ < n_trans; ++trans_){
			     if (clock == "forward"){
			       rates[trans_] = transmod->summary(trans_, s, {tstar}, "hazard")[0];
			     } else { // clock == "mixt"
			       if (transition_types[trans_] == tt_age){
				 rates[trans_] = transmod->summary(trans_, s, {age}, "hazard")[0];
			       } else { // transition_types[trans_] == tt_time (tt_reset is *not* currently available)
				 rates[trans_] = transmod->summary(trans_, s, {tstar}, "hazard")[0];
			       }
			     }
			   } // end loop over transitions
			 } // end loop over states
			 // Rprintf("Utility input calculations\n");
			 // utility input calculations
			 if (valid_R_QALYsStateVal) {
			   int t_index = hesim::hesim_bound(t, obs_index_qalys->time_start_);
			   for (int state_ : live_states) {
			     int obs = (*obs_index_qalys)(k, // strategy
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
			 for (int cost_=0; cost_<n_costs; ++cost_) {
			   int t_index_costs = hesim::hesim_bound(t, obs_index_costs[cost_].time_start_);
			   if (costs_lookup[cost_].method_ == "wlos") {
			     for (int state_ : live_states) {
			       int obs_costs = obs_index_costs[cost_](k, // strategy
								      i, // patient
								      state_,
								      t_index_costs);
			       costs_by_state(state_, cost_) += costs_lookup[cost_].sim(s, obs_costs, type);
			     } // end loop over live states
			   } else if (costs_lookup[cost_].method_ == "starting") {
			     for (int trans_=0; trans_<n_trans; ++trans_) {
			       if (std::count(live_states.begin(), live_states.end(), to[trans_])) {
				 int obs_costs = obs_index_costs[cost_](k, // strategy
									i, // patient
									to(trans_),
									t_index_costs);
				 costs_by_transition(trans_, cost_) = costs_lookup[cost_].sim(s, obs_costs, type);
			       }
			     } // end loop over transitions
			   } else { // costs_lookup[cost_].method_ == "transition"
			     for (int trans_=0; trans_<n_trans; ++trans_) {
			       int obs_costs = obs_index_costs[cost_](k, // strategy
								      i, // patient
								      trans_, // assumes transition
								      t_index_costs);
			       costs_by_transition(trans_, cost_) = costs_lookup[cost_].sim(s, obs_costs, type); // by trans_
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
			 for (int cost_=0; cost_<n_costs; ++cost_) {
			   if (costs_lookup[cost_].method_ == "wlos")
			     dYdt(arma::span((3+cost_)*n_states,(4+cost_)*n_states-1)) += costs_by_state.col(cost_) % Y(arma::span(0,n_states-1))*drr_costs;
			   else if (costs_lookup[cost_].method_ == "starting") {
			     for (int trans_=0; trans_<n_trans; ++trans_)
			       dYdt(to[trans_] + (3+cost_)*n_states) += costs_by_transition(trans_,cost_) * Y(from[trans_]) * rates(trans_) *drr_costs;
			   } else { // costs_lookup[cost_].method_ == "transition"
			     for (int trans_=0; trans_<n_trans; ++trans_)
			       dYdt(from[trans_] + (3+cost_)*n_states) += // costs arbitrarily assigned to the state of origin (from state)
				 costs_by_transition(trans_,cost_) * Y(from(trans_)) * rates(trans_) * drr_costs;
			   }
			 }
		       },
					  p0,
					  times[j-1],
					  times[j],
					  times[j]-times[j-1]);
	    report.row(j) = p0.t();
	  } // end j times_ loop
	  for (int h = 0; h < n_states; ++h){
	    for (int ti = 0; ti < n_times; ++ti){
	      CTSTM_OUT(out,counter);
	      out.t_[counter] = times[ti];
	      out.prob_[counter] = report(ti, h);
	      ++counter;                   
	    } // end cycles loop
	    // report for life-years
	    CTSTM_OUT(out2,counter2);
	    out2.dr_[counter2] = 0.0; // assume no discounting for life-years
	    out2.outcome_[counter2] = "ly";
	    out2.value_[counter2] = report(n_times-1,n_states+h);
	    ++counter2;
	    // report for qalys
	    if (valid_R_QALYsStateVal) {
	      CTSTM_OUT(out2,counter2);
	      out2.dr_[counter2] = dr_qalys;
	      out2.outcome_[counter2] = "qaly";
	      out2.value_[counter2] = report(n_times-1,2*n_states+h);
	      ++counter2;
	    }
	    // report for costs
	    for (int cost_=0; cost_<n_costs; ++cost_) {
	      CTSTM_OUT(out2,counter2);
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

#undef CTSTM_OUT

// test example that does not use hesim
class TestODE {
public:
  typedef arma::vec state_type;
  size_t n_states;
  double discount_rate;
  TestODE(double discount_rate = 0.0)
    : n_states(4), discount_rate(discount_rate) { }
  inline double hweibull(const double t, const double shape, const double scale) {
    return R::dweibull(std::max(1e-100,t),shape,scale,0)/R::pweibull(t,shape,scale,0,0);
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
    for (int j=0; j<4; ++j)
      dYdt(from[j] + 4*n_states) += costs_per_transition[j] * Y(from[j]) * rates(j) * dr;
    // discounted starting costs
    for (int j=0; j<4; ++j)
      dYdt(to[j] + 5*n_states) += starting_costs[to[j]] * Y(from[j]) * rates(j) * dr;
  }
  arma::mat run(arma::vec p0, arma::vec times) {
    using namespace boost::numeric::odeint;
    size_t n_times = times.size();
    // combine the results
    arma::mat combined(n_times,p0.size());
    combined.row(0) = p0.t();
    auto stepper = make_dense_output(1.0e-10, 1.0e-10, runge_kutta_dopri5<state_type>());
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
    return combined;
  }
};
// [[Rcpp::export]]
Rcpp::List runTestODE(arma::vec p0, arma::vec times, double discount_rate = 0.0) {
  using namespace Rcpp;
  return List::create(_("times")=times,
		      _("Y")=TestODE(discount_rate).run(p0,times));
}
