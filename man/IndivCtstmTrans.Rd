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
\item{\code{input_mats}}{Input matrices used to simulate health state transitions by treatment strategy and patient. Must be 
  an object of class \code{\link{input_mats}}. If \code{params} is a list of models, then \code{input_mats} must contain a unique row for each treatment strategy and patient; if \code{params} is a joint model, then \code{input_mats} must contain a unique row for each treatment strategy, patient, and transition.}

\item{\code{params}}{An object of class \code{\link{params_surv}} or \code{\link{params_surv_list}}.}

\item{\code{trans_mat}}{A transition matrix describing the states and transitions in a multi-state model in the format from the \link[mstate]{mstate} package. See the documentation for the argument "trans" in \link[mstate]{msprep}.}

\item{\code{start_state}}{A scalar or vector denoting the starting health state. Default is the first health state. If a vector, must be equal to the number of simulated patients.}

\item{\code{start_time}}{A scalar or vector denoting the time at which to begin the simulation. Default is time 0. If a vector, must be equal to the number of simulated patients.}

\item{\code{clock}}{"reset" for a clock-reset model, "forward" for a clock-forward model, and "mix" for a mixture
of clock-reset and clock-forward models. A clock-reset model is a semi-Markov model in which transition rates depend on time since entering a state. A clock-forward model is a Markov model in which transition rates depend on time since \code{start_time}. If \code{"mix"} is used, then \code{reset_states} must be specified.}

\item{\code{reset_states}}{A vector denoting the states in which time resets. Hazard functions are always a function of elapsed time since either the start of the model or from when time was previously reset. Only used if \code{clock = "mix"}.}

\item{\code{start_age}}{A scalar or vector denoting the starting age of each patient in the simulation. Default is 38. If a vector, must be equal to the number of simulated patients.}

\item{\code{death_state}}{The death state in \code{trans_mat}. Used with \code{max_age} in \code{sim_disease} as patients transition to this state upon reaching maximum age. By default, it is set to the final absorbing state (i.e., a row in \code{trans_mat} with all NAs).}

}}

\section{Methods}{

\describe{

\item{\code{new(input_mats, params, trans_mat, start_state = 1, start_time = 0,
                start_age = 38, death_state = NULL,
                clock = c("reset", "forward", "mix"),
                reset_states = NULL))}}{

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
  
  \item{\code{sim_stateprobs(t, ...)}}{
    Simulate health state probabilities at distinct times.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     \item \code{...}: Additional arguments to pass to \code{sim_disease}.
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

\seealso{
\code{\link{create_IndivCtstmTrans}}, \code{\link{IndivCtstm}}
}

\examples{
library("flexsurv")

# Simulation data
strategies <- data.frame(strategy_id = c(1, 2, 3))
patients <- data.frame(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))

# Multi-state model with transition specific models
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))
fits <- vector(length = max(tmat, na.rm = TRUE), mode = "list")
for (i in 1:length(fits)){
  fits[[i]] <- flexsurvreg(Surv(years, status) ~ 1, 
                            data = bosms3[bosms3$trans == i, ],
                           dist = "exp")
}
fits <- flexsurvreg_list(fits)

# Simulation model
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)
fits_data <- expand(hesim_dat)
transmod <- create_IndivCtstmTrans(fits, input_data = fits_data, 
                                    trans_mat = tmat,
                                    n = 2,
                                    point_estimate = FALSE)       
head(transmod$hazard(c(1, 2, 3)))
head(transmod$cumhazard(c(1, 2, 3)))
transmod$sim_stateprobs(t = c(0, 5, 10))[t == 5]
}