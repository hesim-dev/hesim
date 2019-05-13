\docType{data}
\name{StateVals}
\alias{StateVals}

\title{Model for state values}
\format{\code{\link{R6Class}} object.}
\usage{
StateVals
}
\description{
Simulate values (i.e., utility or costs) associated with health states in a partitioned survival or state transition model.
}

\section{Fields}{

  \describe{

  \item{\code{input_mats}}{Input matrices used to simulate state values. Must be an object of class \code{\link{input_mats}} where each row in the matrix "X" is a unique strategy, patient, health state, and (optionally) time interval. "X" must be sorted by strategy, patient, health state, and (if specified) time interval.}

  \item{\code{params}}{An object of class \code{\link{params_lm}} or \code{\link{params_mean}}, which contain the parameters for simulatingstate values.}
  
  \item{\code{max_t}}{A vector specifying the maximum duration that state values are computed for in the
simulation once a patient has entered a (new) health state. For example, if a drug was only given for 1 year, then \code{max_t} would be set to 1 (year) and no costs would accrue after that year. The length should be equal to the number of rows in \code{input_mats}. If not specified, then \code{max_t} is implictly assumed to equal infinity. Can only be used in individual-level simulations (i.e., \code{\link{IndivCtstm}}), which are able to track time spent in each state.}

  }
}

\section{Methods}{

\describe{

\item{\code{new(input_mats, params)}}{

  Constructor for the class.

  }

    \item{\code{sim(t, type = c("predict", "random"))}}{

    Simulate state values with either predicted means or random samples by treatment strategy, patient, health state, and (optionally) time interval.
    
    \itemize{
     \item \code{t}: A numeric vector of times. The first element must be 0.
     \item \code{type}: \code{"predict"} for mean values or \code{"random"} for random samples.
     }

    }

    \item{\code{check()}}{

    Input validation for class. Checks that fields are the correct type. 

    }

  }

}

\examples{
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

# Create StateVals object
fit_costs_medical <- stats::lm(costs ~ female + state_name, data = psm4_exdata$costs$medical)
dat <- expand(hesim_dat, by = c("strategies", "patients", "states"))
costs_medical <- create_StateVals(fit_costs_medical, input_data = dat, n = 5)

# Predict
head(costs_medical$sim(t = 5, type = "predict"))
}