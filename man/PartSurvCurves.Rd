\docType{data}
\name{PartSurvCurves}
\alias{PartSurvCurves}

\title{Partitioned survival curves}
\format{\code{\link{R6Class}} object.}
\usage{
PartSurvCurves
}
\description{
Summarize \code{n-1} survival curves for an \code{n} state partitioned survival model.

}
\section{Fields}{

\describe{
\item{\code{data}}{Data used to predict state values by strategy and patient. Must be 
  an object of class \code{\link{input_data}} where each row of a matrix in "X" is a unique strategy and patient. All   matrices in "X" must be sorted by strategy and patient.}

\item{\code{params}}{An object of class \code{\link{params_surv_list}}.}

}}

\section{Methods}{

\describe{

\item{\code{new(data, params)}}{

  Constructor for the class.

  \itemize{
    \item{\code{data}}{Data used to summarize survival curves by strategy and patient. Must be 
  an object of class \code{\link{input_data}} where each row in "X" is a strategy and patient pair.}

  \item{\code{params}}{An object of class \code{\link{partsurvfit}}, which contains the parameters for predicting
  survival.}
  }
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

    Input validation for class; checks that fields are the correct type. 

    }

}


}

\examples{
library("flexsurv")

# Simulation data
dt.strategies <- data.frame(strategy_id = c(1, 2, 3))
dt.patients <- data.frame(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
hesim.dat <- hesim_data(strategies = dt.strategies,
                              patients = dt.patients)

# Fit survival models
surv.data <- part_surv4_simdata$survival
fit1 <- flexsurv::flexsurvreg(Surv(endpoint1_time, endpoint1_status) ~ age,
                              data = surv.data, dist = "exp")
fit2 <- flexsurv::flexsurvreg(Surv(endpoint2_time, endpoint2_status) ~ age,
                              data = surv.data, dist = "exp")
fit3 <- flexsurv::flexsurvreg(Surv(endpoint3_time, endpoint3_status) ~ age,
                              data = surv.data, dist = "exp")
fits <- partsurvfit(flexsurvreg_list(fit1, fit2, fit3), data = surv.data)

# Form PartSurvCurves
curves.data <- expand_hesim_data(hesim.dat, by = c("strategies", "patients"))
part.surv.curves <- form_PartSurvCurves(fits, data = curves.data, n = 3)

# Summarize survival curves
quantiles <- part.surv.curves$quantile(p = c(.25, .5, .75))
head(quantiles)

survival <- part.surv.curves$survival(t = seq(0, 3, by = .1))
head(survival)

rmst <- part.surv.curves$rmst(t = c(2, 5))
head(rmst)
}
