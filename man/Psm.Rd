\docType{data}
\name{Psm}
\alias{Psm}

\title{N-state partitioned survival model}
\format{\code{\link{R6Class}} object.}
\usage{
Psm
}
\description{
Simulate outcomes from an N-state partitioned survival model.

}
\section{Fields}{

\describe{
\item{\code{survival_models}}{The survival models used to predict survival curves. Must be
an object of class \code{\link{PsmCurves}}.}

\item{\code{utility_model}}{The model used to predict utility by health state. Must be
an object of class \code{\link{StateVals}}.}

\item{\code{cost_models}}{The models used to predict costs by health state. Must be a list of
objects of class \code{\link{StateVals}}, where each element of the list represents a 
different cost category.}

\item{\code{n_states}}{Number of states in the partitioned survival model.}

\item{\code{t_}}{A numeric vector of times at which survival curves were predicted. Determined
by the argument \code{t} in \code{sim_curves}.}

\item{\code{survival_}}{Survival curves generated using \code{sim_curves}.}

\item{\code{stateprobs_}}{Health state probabilities as a function of time generated
using \code{sim_stateprobs}.}

\item{\code{costs_}}{Total (discounted) costs by type generated using \code{sim_costs}. }

\item{\code{qalys_}}{Total (discounted) quality-adjusted life-years (QALYs) generated using \code{sim_qalys}.}

}}

\section{Methods}{

\describe{

\item{\code{new(survival_models, utility_model = NULL, cost_models = NULL)}}{

  Constructor for the class. Note that the number of health states, \code{n_states}, is set equal to the number of survival
  models plus one. 
  }
  
  \item{\code{sim_survival(t)}}{
    Simulate survival curves as a function of time. Equivalent to the member function
    \code{survival} in \code{\link{PsmCurves}}. 
    
     \itemize{
     \item \code{t}: A numeric vector of times. The first element must be 0.
     }
  }
  
  \item{\code{sim_stateprobs()}}{
    Simulate the probability of being in each of N health states using the survival curves generated 
    from \code{sim_curves}. 
  }
  
    \item{\code{sim_qalys(dr = .03)}}{
    Simulate (discounted) QALYs over the times selected in \code{t} associated with each health state based on the state probabilities calculated using \code{sim_stateprobs}. See "Details".
    
      \itemize{
     \item \code{dr}: Discount rate to apply to QALYs. May be a vector in which case QALYs are calculated for each element in \code{dr}.
     }
  }

  
   \item{\code{sim_costs(dr = .03)}}{
    Simulate (discounted) costs for each cost type over the times selected in \code{t} associated with each health state based on the state probabilities calculated using \code{sim_stateprobs}. See "Details".
    
      \itemize{
     \item \code{dr}: Discount rate to apply to costs. May be a vector in which case costs for each cost type are     calculated for each element in \code{dr}.
       }
  }

  \item{\code{check()}}{
 
    Input validation for class. Checks that fields are the correct type. 
  }

}
}

\section{Details}{
Discounted costs and QALYs are calculated by integrating the "weighted" probability of being in each state. 
Weights are a function of the discount factor and the state value predicted using either the cost or QALY model. 
Mathematically, discounted costs and QALYs in health state \eqn{s} are calculated as,

\deqn{\int_0^T w_h e^{-rt} P_h(t) dt },

where for health state \eqn{h} and time {t}, \eqn{w_h} is the predicted cost or QALY weight, \eqn{r} is the discount rate, and \eqn{P_h(t)} is the probability of being in a given health state. The integral is calculated numerically using the composite trapezoid rule from the points in \code{t_}.

}


\examples{
library("flexsurv")

# Simulation data
dt_strategies <- data.frame(strategy_id = c(1, 2, 3))
dt_patients <- data.frame(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
dt_states <- data.frame(state_id =  seq(1, 3),
                        state_name = paste0("state", seq(1, 3)))
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients,
                        states = dt_states)
N <- 3

# Survival models
surv_data <- psm4_exdata$survival
fit1 <- flexsurv::flexsurvreg(Surv(endpoint1_time, endpoint1_status) ~ age,
                              data = surv_data, dist = "exp")
fit2 <- flexsurv::flexsurvreg(Surv(endpoint2_time, endpoint2_status) ~ age,
                              data = surv_data, dist = "exp")
fit3 <- flexsurv::flexsurvreg(Surv(endpoint3_time, endpoint3_status) ~ age,
                              data = surv_data, dist = "exp")
fits <- partsurvfit(flexsurvreg_list(fit1, fit2, fit3), data = surv_data)

curves_data <- expand_hesim_data(hesim_dat, by = c("strategies", "patients"))
psm_curves <- create_PsmCurves(fits, data = curves_data, n = N)

# Cost models
expanded_data <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "states"))
fit_costs_medical <- stats::lm(costs ~ female + state_name, data = psm4_exdata$costs$medical)
psm_costs_medical <- create_StateVals(fit_costs_medical, data = expanded_data, n = N)

# Utility model
psm_utility_data <- create_input_data(formula_list(mu = ~1), expanded_data)
psm_utility <- StateVals$new(data = psm_utility_data,
                             params = params_lm(coef = runif(N, .6, .8)))

# Partitioned survival decision model
psm <- Psm$new(survival_models = psm_curves,
               utility_model = psm_utility,
               cost_models = list(medical = psm_costs_medical))
times <- seq(0, 5, .05)
psm$sim_survival(t = seq(0, 5, .05))
psm$sim_stateprobs()
psm$sim_costs(dr = .03)
head(psm$costs_)
head(psm$sim_qalys(dr = .03)$qalys_)
}
