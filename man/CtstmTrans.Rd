\docType{data}
\name{CtstmTrans}
\alias{CtstmTrans}

\title{Health state transitions for continuous time state transition models}
\format{\code{\link{R6Class}} object.}
\usage{
CtstmTrans
}
\description{
Summarize health state transitions in a continuous time state transition model with parameters that were estimated using a multi-state model.

}
\section{Fields}{

\describe{
\item{\code{data}}{Data used to predict health state transitions by strategy and patient. Must be 
  an object of class \code{\link{input_data}}.}

\item{\code{params}}{An object of class \code{\link{params_surv_list}}.}

\item{\code{trans_mat}}{A transition matrix describing the states and transitions in a multi-state model in the format from the \link[mstate]{mstate} package. See the documentation for the argument "trans" in \link[mstate]{msprep}.}

}}

\section{Methods}{

\describe{

\item{\code{new(data, params)}}{

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
  }

    
    \item{\code{check()}}{

    Input validation for class. Checks that fields are the correct type. 

    }

}


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
}