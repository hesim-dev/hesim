\docType{data}
\name{IndivCtstm}
\alias{IndivCtstm}

\title{Individual-level continuous time state transition model}
\format{\code{\link{R6Class}} object.}
\usage{
IndivCtstm
}
\description{
Simulate outcomes from an individual-level continuous time state transition model (CTSTM). The class is currently only valid for semi-Markov (i.e., "clock-reset") models; however,
functionality will soon be added for Markov  (i.e., "clock-forward") models and mixtures of clock-forward and clock-reset models.

}
\section{Fields}{

\describe{
\item{\code{trans_model}}{Model for health state transitions. Must be an object of class \code{\link{IndivCtstmTrans}}.}

\item{\code{utility_model}}{The model used to predict utility by health state. Must be an object of class \code{\link{StateVals}}.}

\item{\code{cost_models}}{The models used to predict costs by health state. Must be a list of objects of class \code{\link{StateVals}}, where each element of the list represents a different cost category.}

\item{\code{disprog_}}{A \code{\link[data.table]{data.table}} simulated using \code{sim_disease} containing the following columns:

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

}

\item{\code{stateprobs_}}{A \code{\link[data.table]{data.table}} of health state probabilities as a function of time simulated using \code{sim_stateprobs}. See the description of the output of \code{sim_stateprobs} in \code{\link{IndivCtstmTrans}}.

}

\item{\code{qalys_}}{A \code{\link[data.table]{data.table}} of quality-adjusted life-years (QALYs) simulated using \code{sim_qalys}. Columns are:

\describe{
\item{sample}{A random sample from the probabilistic sensitivity analysis (PSA).}
\item{strategy_id}{The treatment strategy ID.}
\item{patient_id}{A patient ID.}
\item{dr}{The discount rate.}
\item{qalys}{Simulated QALYs.}
\item{lys}{Simulated life-years.}
}

Note that the \code{lys} column is only computed if \code{lys = TRUE}.

}

\item{\code{costs_}}{A \code{\link[data.table]{data.table}} of costs by category simulated using \code{sim_costs}. Columns are:

\describe{
\item{sample}{A random sample from the PSA.}
\item{strategy_id}{The treatment strategy ID.}
\item{patient_id}{A patient ID.}
\item{dr}{The discount rate.}
\item{category}{The cost category.}
\item{costs}{Simulated costs.}
}

}

}}

\section{Methods}{

\describe{

\item{\code{new(trans_model = NULL, utility_model = NULL, cost_models = NULL)}}{

  Constructor for the class.
  
}
  
  \item{\code{sim_disease(max_t = 100, max_age = 100)}}{
    Simulate disease progression.
    
   \itemize{
     \item \code{max_t}: A scalar or vector denoting the number of time periods to simulate the model. If a vector, must be equal to the number of simulated patients. 
    \item \code{max_age}: A scalar or vector denoting the maximum age to simulate each patient until. If a vector, must be equal to the number of simulated patients.
   }
     
    Returns an instance of self with simulated output stored in \code{disprog_}. 
  }
  
  \item{\code{sim_stateprobs(t)}}{
    Simulate the probability of being in each health state as a function of time using the simulation output stored in \code{disprog_}.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     }
     
    Returns an instance of self with simulated output stored in \code{stateprobs_}.  
  }
  
  \item{\code{sim_qalys(dr = .03, type = c("predict", "random"), by_patient = FALSE, lys = TRUE)}}{
    Compute simulated (mean discounted) QALYs using the simulation output stored in \code{disprog_} by random sample from the PSA, treatment strategy, health state, and (optionally) patient.
    
     \itemize{
     \item \code{dr}: Discount rate to apply to QALYs. May be a vector in which case QALYs are calculated for each element in \code{dr}.
     \item \code{type}: \code{predict} for mean values or \code{random} for random samples as in \code{$sim()} in \code{\link{StateVals}}.
     \item \code{by_patient}: If \code{TRUE}, QALYs are computed at the patient level. If \code{FALSE}, QALYs are averaged across patients by
     health state.
     \item \code{lys}: If \code{TRUE}, then life-years are computed in addition to QALYs.
     }
     
     Returns an instance of self with simulated output stored in \code{qalys_}.
  } 
  
  \item{\code{sim_costs(dr = .03, type = c("predict", "random"), by_patient = FALSE, max_t = Inf)}}{
    Compute simulated (mean discounted) costs using the simulation output stored in \code{disprog_} by random sample from the PSA, treatment strategy, health state, and 
    (optionally) patient.
    
     \itemize{
     \item \code{dr}: Discount rate to apply to costs. May be a vector in which case costs are calculated for each element in \code{dr}.
      \item \code{type}: \code{predict} for mean values or \code{random} for random samples as in \code{$sim()} in \code{\link{StateVals}}.
     \item \code{by_patient}: If \code{TRUE}, costs are computed at the patient level. If \code{FALSE}, costs are averaged across patients by
     health state.      
      \item \code{max_t}: Maximum time duration to compute costs once a patient has entered a (new) health state. By default, equal to \code{Inf}, 
      so that costs are computed over the entire duration that a patient is in a given health state. If time varies by each cost category, then time can also be passed
      as a numeric vector of length equal to the number of cost categories (e.g., \code{c(1, 2, Inf, 3)} for a model with 4 cost categories).
     }
     
    Returns an instance of self with simulated output stored in \code{costs_}. 
  },
  
  \item{\code{summarize()}}{
    Produce a summary measure of costs and QALYs by random sample from the PSA and treatment strategy.
     
     Returns an object of summary cost-effectiveness output of class \code{\link{ce}}. 
  }    
  

}

}

\seealso{
\code{\link{create_IndivCtstmTrans}}, \code{\link{IndivCtstmTrans}}
}


\examples{
library("flexsurv")

# Treatment strategies, target population, and model structure
strategies <- data.frame(strategy_id = c(1, 2))
patients <- data.frame(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
states <- data.frame(state_id = c(1, 2))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)

# Parameter estimation
## Multi-state model
tmat <- rbind(c(NA, 1, 2),
              c(3, NA, 4),
              c(NA, NA, NA))
fits <- vector(length = max(tmat, na.rm = TRUE), mode = "list")
surv_dat <- data.frame(ctstm3_exdata$transitions)
for (i in 1:length(fits)){
  fits[[i]] <- flexsurvreg(Surv(years, status) ~ factor(strategy_id), 
                           data = surv_dat,
                           subset = (trans == i),
                           dist = "weibull")
}
fits <- flexsurvreg_list(fits)

## Utility
utility_tbl <- stateval_tbl(data.frame(state_id = states$state_id,
                                       mean = ctstm3_exdata$utility$mean,
                                       se = ctstm3_exdata$utility$se),
                            dist = "beta",
                            hesim_data = hesim_dat)
## Costs
drugcost_tbl <- stateval_tbl(data.frame(strategy_id = strategies$strategy_id,
                                       est = ctstm3_exdata$costs$drugs$costs),
                            dist = "fixed",
                            hesim_data = hesim_dat) 
medcost_tbl <- stateval_tbl(data.frame(state_id = states$state_id,
                                       mean = ctstm3_exdata$costs$medical$mean,
                                       se = ctstm3_exdata$costs$medical$se),
                            dist = "gamma",
                            hesim_data = hesim_dat)  

# Economic model
n_samples = 2

## Construct model
### Transitions 
transmod_data <- expand(hesim_dat)
transmod <- create_IndivCtstmTrans(fits, data = transmod_data, trans_mat = tmat,
                                   n = n_samples)
                        
### Utility 
utilitymod <- create_StateVals(utility_tbl, n = n_samples)

### Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples)
costmods <- list(drugs = drugcostmod,
                 medical = medcostmod)
                 
### Combine
ictstm <- IndivCtstm$new(trans_model = transmod,
                         utility_model = utilitymod,
                         cost_models = costmods)


## Simulate outcomes
head(ictstm$sim_disease()$disprog_)
head(ictstm$sim_stateprobs(t = c(0, 5, 10))$stateprobs_[t == 5])
ictstm$sim_qalys(dr = .03)
ictstm$sim_costs(dr = .03)
head(ictstm$summarize())
}