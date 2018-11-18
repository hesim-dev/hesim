\docType{data}
\name{PsmCurves}
\alias{PsmCurves}

\title{Partitioned survival curves}
\format{\code{\link{R6Class}} object.}
\usage{
PsmCurves
}
\description{
Summarize \code{n-1} survival curves for an \code{N} state partitioned survival model.

}
\section{Fields}{

\describe{
\item{\code{input_mats}}{Input matrices used to predict state values by strategy and patient. Must be 
  an object of class \code{\link{input_mats}} where each row of a matrix in "X" is a unique strategy and patient. All   matrices in "X" must be sorted by strategy and patient.}

\item{\code{params}}{An object of class \code{\link{params_surv_list}}.}

}}

\section{Methods}{

\describe{

\item{\code{new(input_mats, params)}}{

  Constructor for the class.
  
}
  
  \item{\code{hazard(t)}}{
    Predict the hazard function as a function of time.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     }
  }
  
  \item{\code{cumhazard(t)}}{
    Predict the cumulative hazard function as a function of time.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     }
  }
  
  \item{\code{survival(t)}}{
    Predict the survival function as a function of time.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     }
  }
  
  \item{\code{rmst(t, dr = 0)}}{
    Predict (discounted) restricted mean survival time.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     \item \code{dr}: Discount rate. 
     }
  }
  
    \item{\code{quantile(p)}}{
    Predict quantiles of the survival distributions. 
    
     \itemize{
     \item \code{p}: A numeric vector of probabilities for calculating quantiles.
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
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients)

# Fit survival models
surv_est_data <- psm4_exdata$survival
fit1 <- flexsurv::flexsurvreg(Surv(endpoint1_time, endpoint1_status) ~ age,
                              data = surv_est_data, dist = "exp")
fit2 <- flexsurv::flexsurvreg(Surv(endpoint2_time, endpoint2_status) ~ age,
                              data = surv_est_data, dist = "exp")
fit3 <- flexsurv::flexsurvreg(Surv(endpoint3_time, endpoint3_status) ~ age,
                              data = surv_est_data, dist = "exp")
fits <- flexsurvreg_list(fit1, fit2, fit3)

# Form PsmCurves
surv_input_data <- expand(hesim_dat, by = c("strategies", "patients"))
psm_curves <- create_PsmCurves(fits, data = surv_input_data, n = 3,
                                bootstrap = TRUE, est_data = surv_est_data)

# Summarize survival curves
head(psm_curves$quantile(p = c(.25, .5, .75)))
head(psm_curves$survival(t = seq(0, 3, by = .1)))
head(psm_curves$rmst(t = c(2, 5)))

}
