\docType{data}
\name{StateVals}
\alias{StateVals}

\title{State values}
\format{\code{\link{R6Class}} object.}
\usage{
StateVals
}
\description{
Simulate values (i.e., utility or costs) associated with health states in a partitioned survival or state transition model.
}

\section{Fields}{

  \describe{

  \item{\code{data}}{Data used to simulate state values by strategy, patient, and health state. Must be 
  an object of class \code{\link{input_data}} where each row in the matrix "X" is a unique strategy, patient, and
  health state. "X" must be sorted by strategy, patient, and health state.}

  \item{\code{params}}{An object of class \code{\link{params_lm}}, which contains the parameters for simulating
  state values.}

  }
}

\section{Methods}{

\describe{

\item{\code{new(data, params)}}{

  Constructor for the class.

  }

    \item{\code{sim(t, type = c("predict", "random"))}}{

    Simulate state values with either predicted means or random samples by treatment strategy, patient, and health state.
    
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
dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "states"))
costs_medical <- create_StateVals(fit_costs_medical, data = dat, n = 5)

# Predict
head(costs_medical$sim(t = 5, type = "predict"))
}