\docType{data}
\name{IndivCtstm}
\alias{IndivCtstm}

\title{Individual-level continuous time state transition model}
\format{\code{\link{R6Class}} object.}
\usage{
IndivCtstm
}
\description{
Simulate outcomes from an individual-level continuous time state transition model (CTSTM). The class is only valid for semi-Markov (i.e., "clock-reset") models.

}
\section{Fields}{

\describe{
\item{\code{trans_model}}{Model for health state transitions. Must be an object of class \code{\link{CtstmTrans}}.}

\item{\code{utility_model}}{The model used to predict utility by health state. Must be an object of class \code{\link{StateVals}}.}

\item{\code{cost_models}}{The models used to predict costs by health state. Must be a list of objects of class \code{\link{StateVals}}, where each element of the list represents a different cost category.}

\item{\code{disease_prog_}}{Disease progression simulated using \code{sim_disease}; that is, a tidy  \code{\link[data.table]{data.table}} containing simulated trajectories through the multi-state model. Columns include: 

\describe{
\item{sample}{A random sample from the PSA.}
\item{strategy_id}{The treatment strategy ID.}
\item{patient_id}{The patient ID.}
\item{from}{The health state ID transitioned from.}
\item{to}{The health state ID transitioned to.}
\item{final}{An indicator equal to 1 if a patient is in their final health state during the simulation and 0 otherwise.}
\item{time_start}{The time at the start of the interval.}
\item{time_start}{The time at the end of the interval.}
}
}

\item{\code{stateprobs_}}{A tidy \code{\link[data.table]{data.table}} of health state probabilities as a function of time simulated using \code{sim_stateprobs}. Columns include:

\describe{
\item{sample}{A random sample from the PSA.}
\item{strategy_id}{The treatment strategy ID.}
\item{state_id}{The health state ID.}
\item{t}{The time at which a state probability is computed.}
\item{prob}{The probability of being in a given health state.}
}

}

}}

\section{Methods}{

\describe{

\item{\code{new(trans_model, utility_model = NULL, cost_models = NULL)}}{

  Constructor for the class.
  
}
  
  \item{\code{sim_disease(max_t = 100, max_age = 100)}}{
    Simulate disease progression.
    
     \itemize{
     \item \code{max_t}: Maximum time to simulate disease until.
     \item \code{max_age}: Maximum age that a patient can live. Upon reaching this age the simulation stops and a patient moves to the "death state" specified by \code{death_state}.
     }
     
    The simulated output is stored in \code{disease_prog_}. 
  }
  
  \item{\code{sim_stateprobs(t)}}{
    Simulate the probability of being in each health state as a function of time.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     }
     
    The simulated output is stored in \code{stateprobs_}.  
  }

}

}

\seealso{
\code{\link{create_CtstmTrans}}, \code{\link{CtstmTrans}}
}


\examples{
library("flexsurv")
library("data.table")

# Simulation data
dt_strategies <- data.frame(strategy_id = c(1, 2, 3))
dt_patients <- data.frame(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))

# Multi-state model
fits <- vector(length = 3, mode = "list")
dat <- data.table(bosms3)
for (i in 1:length(fits)){
  fits[[i]] <- flexsurvreg(Surv(years, status) ~ 1, data = dat[trans == i],
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
ictstm <- IndivCtstm$new(trans_model = transmod)
head(ictstm$sim_disease()$disease_prog_)
ictstm$sim_stateprobs(t = c(0, 5, 10))$stateprobs_[t == 5]
}