\docType{data}
\name{IndivCtstmTrans}
\alias{IndivCtstmTrans}

\title{Transitions for an individual-level continuous time state transition model}
\format{\code{\link{R6Class}} object.}
\usage{
IndivCtstmTrans
}
\description{
Summarize health state transitions in an individual-level continuous time state transition model with parameters that were estimated using a multi-state model.

}
\section{Fields}{

\describe{
\item{\code{data}}{Data used to simulate health state transitions by treatment strategy and patient. Must be 
  an object of class \code{\link{input_data}}. If \code{params} is a list of models, then \code{data} must contain a unique row for each treatment strategy and patient; if \code{params} is a joint model, then \code{data} must contain a unique row for each treatment strategy, patient, and transition.}

\item{\code{params}}{An object of class \code{\link{params_surv}} or \code{\link{params_surv_list}}.}

\item{\code{trans_mat}}{A transition matrix describing the states and transitions in a multi-state model in the format from the \link[mstate]{mstate} package. See the documentation for the argument "trans" in \link[mstate]{msprep}.}

\item{\code{start_state}}{The vectorized (see "Details") starting health state. Default is the first health state.}

\item{\code{start_time}}{The vectorized (see "Details")  time at which to begin the simulation. Default is time 0.}

\item{\code{start_age}}{The vectorized (see "Details") starting age of each patient in the simulation. Default is 38.}

\item{\code{death_state}}{The death state in \code{trans_mat}. Used with \code{max_age} in \code{sim_disease} as patients transition to this state upon reaching maximum age. By default, it is set to the final absorbing state (i.e., a row in \code{trans_mat} with all NAs).}

}}

\section{Methods}{

\describe{

\item{\code{new(data, params, trans_mat, start_state = 1, start_time = 0,
                start_age = 38, death_state = NULL)}}{

  Constructor for the class.
  
}
  
  \item{\code{hazard(t)}}{
    Predict the hazard functions for each health state transition.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     }
  }
  
  \item{\code{cumhazard(t)}}{
    Predict the cumulative hazard functions for each health state transition.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     }
  },
  
  \item{\code{sim_disease(max_t = 100, max_age = 100)}}{
  Simulate disease progression.
  
   \itemize{
     \item \code{max_t}: The vectorized (see "Details") maximum number of time periods to simulate the model.
    \item \code{max_age}: Maximum age to simulate each patient until. Default is 100.
   }
   
   Returns an object of class "indiv_ctstm_disprog". See "Details". 
  },
  
  
  \item{\code{sim_stateprobs(disprog = NULL, t, ...)}}{
    Simulate health state probabilities at distinct times.
    
     \itemize{
     \item \code{disprog}: If \code{NULL}, then state probabilities are simulated by first simulating disease progression using \code{sim_disease}. Otherwise,
     state probabilities are computed using the simulation output from \code{disprog}. Must be of class "indiv_ctstm_disprog".
     \item \code{t}: A numeric vector of times.
     \item \code{...}: Additional arguments to pass to \code{sim_disease} when \code{disprog = NULL}.
     }
     
    Returns a \code{\link[data.table]{data.table}} with the following columns:
    
    \describe{
    \item{sample}{A random sample from the PSA.}
    \item{strategy_id}{The treatment strategy ID.}
    \item{state_id}{The health state ID.}
    \item{t}{The time at which a state probability is computed.}
    \item{prob}{The probability of being in a given health state.}
    }
  }  

    
    \item{\code{check()}}{

    Input validation for class. Checks that fields are the correct type. 

    }

}


}

\section{Details}{
An object of class "indiv_ctstm_disprog" is a list containing the following elements:
\describe{
\item{sim}{A \code{\link[data.table]{data.table}} of simulated trajectories through a multi-state model.}
\item{n_samples}{The number of random samples of the parameters used in the probabilistic sensitivity analysis.}
\item{n_states}{The number of health states in the model.}
\item{unique_strategy_id}{The ID's of the simulated treatment strategies.}
\item{unique_patient_id}{The ID's of the simulated patients.}
}

The \code{\link[data.table]{data.table}} \code{sim} contains the following columns:

\describe{
\item{sample}{A random sample from the PSA.}
\item{strategy_id}{The treatment strategy ID.}
\item{patient_id}{The patient ID.}
\item{from}{The health state ID transitioned from.}
\item{to}{The health state ID transitioned to.}
\item{final}{An indicator equal to 1 if a patient is in their final health state during the simulation and 0 otherwise.}
\item{time_start}{The time at the start of the interval.}
\item{time_stop}{The time at the end of the interval.}
}

The fields \code{start_state}, \code{start_time}, \code{start_age}, and the argument \code{max_t} in \code{sim_disease} are vectorized across patients, treatment strategies, and random samples of the parameters. Each argument can be passed in three ways:

\itemize{
\item As a scalar, implying that values are constant across patients, treatment strategies, and samples.
\item As a vector of length equal to the number of patients, implying that values vary across patients but are constant across treatment strategies and samples.
\item As a vector of length equal to the product of the number of patients, treatment strategies, and samples. In this case the order should be sorted by sample, treatment strategy, and patient (as with an object of class "indiv_ctstm_disprog"). This option can be used to combine results from multiple simulations using \code{sim_disease}. Specifically, values previously simulated can be used as arguments to simulations from a new instance of \code{IndivCtstmTrans}. 
}
}

\seealso{
\code{\link{create_IndivCtstmTrans}}, \code{\link{IndivCtstm}}
}

\examples{
library("flexsurv")

# Simulation data
dt_strategies <- data.frame(strategy_id = c(1, 2, 3))
dt_patients <- data.frame(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))

# Multi-state model with transition specific models
fits <- vector(length = 3, mode = "list")
dat <- data.frame(bosms3)
for (i in 1:length(fits)){
  fits[[i]] <- flexsurvreg(Surv(years, status) ~ 1, data = dat[dat$trans == i, ],
                           dist = "exp")
}
fits <- flexsurvreg_list(fits)
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))

# Simulation model
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients)
fits_data <- expand_hesim_data(hesim_dat)
transmod <- create_IndivCtstmTrans(fits, data = fits_data, trans_mat = tmat,
                                    n = 2,
                                    point_estimate = FALSE)       
head(transmod$hazard(c(1, 2, 3)))
head(transmod$cumhazard(c(1, 2, 3)))
disprog <- transmod$sim_disease()
names(disprog)
head(disprog$sim)
transmod$sim_stateprobs(t = c(0, 5, 10))[t == 5]
}