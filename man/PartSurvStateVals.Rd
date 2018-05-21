\docType{data}
\name{PartSurvStateVals}
\alias{PartSurvStateVals}

\title{State values for partitioned survival model}
\format{\code{\link{R6Class}} object.}
\usage{
PartSurvStateVals
}
\description{
Predict or values (i.e., utility or costs) associated with health states in a partitioned survival model.
}

\section{Fields}{

  \describe{

  \item{\code{data}}{Data used to predict state values by strategy, patient, and health state. Must be 
  an object of class \code{\link{input_data}} where each row in the matrix "X" is a unique strategy, patient, and
  health state. "X" must be sorted by strategy, patient, and health state.}

  \item{\code{params}}{An object of class \code{\link{params_lm}}, which contains the parameters for predicting
  state values.}

  }
}

\section{Methods}{

\describe{

\item{\code{new(data, params)}}{

  Constructor for the class.

  }

    \item{\code{predict()}}{

    Predict mean state values by treatment strategy, patient, and health state.

    }

    \item{\code{check()}}{

    Input validation for class. Checks that fields are the correct type. 

    }

  }

}

\examples{
# Simulation data
dt.strategies <- data.frame(strategy_id = c(1, 2, 3))
dt.patients <- data.frame(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
dt.states <- data.frame(state_id =  seq(1, 3),
                           state_name = paste0("state", seq(1, 3)))
hesim.dat <- hesim_data(strategies = dt.strategies,
                              patients = dt.patients,
                              states = dt.states)

# Create PartSurvStateVals object
fit.costs.medical <- stats::lm(costs ~ female + state_name, data = part_surv4_simdata$costs$medical)
dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", "states"))
part.surv.costs.medical <- form_PartSurvStateVals(fit.costs.medical, data = dat, n = 5)

# Predict
head(part.surv.costs.medical$predict())
}