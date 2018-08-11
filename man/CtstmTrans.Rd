\docType{data}
\name{CtstmTrans}
\alias{CtstmTrans}

\title{Transitions for continuous time state transition model}
\format{\code{\link{R6Class}} object.}
\usage{
CtstmTrans
}
\description{
Summarize health state transitions in a continuous time state transition model with parameters that were estimated using a multi-state model.

}
\section{Fields}{

\describe{
\item{\code{data}}{Data used to simulate health state transitions by treatment strategy and patient. Must be 
  an object of class \code{\link{input_data}}. If \code{params} is a list of models, then \code{data} must contain a unique row for each treatment strategy and patient; if \code{params} is a joint model, then \code{data} must contain a unique row for each treatment strategy, patient, and transition.}

\item{\code{params}}{An object of class \code{\link{params_surv}} or \code{\link{params_surv_list}}.}

\item{\code{trans_mat}}{A transition matrix describing the states and transitions in a multi-state model in the format from the \link[mstate]{mstate} package. See the documentation for the argument "trans" in \link[mstate]{msprep}.}

\item{\code{start_ages}}{The starting age of each patient in the simulation. Must be the same length as the number of patients in \code{data}. Used alongside \code{max_age} in \code{sim_stateprobs()} with \code{type = "individual"}; does not otherwise need to be specified.}

\item{\code{death_state}}{The death state in \code{trans_mat}. Used with \code{max_age} in \code{sim_stateprobs} with \code{type = "individual"} and in \code{IndivCtstm$sim_disease} as patients transition to this state upon reaching maximum age; does not otherwise need to be specified. By default, it is set to the final absorbing state (i.e., a row in \code{trans_mat} with all NAs).}

}}

\section{Methods}{

\describe{

\item{\code{new(data, params, trans_mat, start_ages = rep(38, data$n_patients),
                          death_state = NULL)}}{

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
  
  \item{\code{sim_stateprobs(t, start_state = 1, max_t = 100, max_age = 100,
                            clock = "reset")}}{
    Simulate health state probabilities as a function of time.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     \item \code{start_state}: The starting health state. Default is the first
     health state. 
     \item \code{max_t}: Maximum number of time periods to simulate the model for. Default is 100.
     \item \code{max_age}: Maximum age to simulate each patient until. Default is 100.
     \item \code{clock}: Currently only supports \code{"reset"}, which implies that the time-scale of the
     cumulative hazards depends on time in the present state. 
     }
  }  

    
    \item{\code{check()}}{

    Input validation for class. Checks that fields are the correct type. 

    }

}


}

\seealso{
\code{\link{create_CtstmTrans}}, \code{\link{IndivCtstm}}
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
transmod <- create_CtstmTrans(fits, data = fits_data, trans_mat = tmat,
                              point_estimate = TRUE)       
head(transmod$hazard(c(1, 2, 3)))
head(transmod$cumhazard(c(1, 2, 3)))
head(transmod$sim_stateprobs(t = c(0, 5, 10)))
}