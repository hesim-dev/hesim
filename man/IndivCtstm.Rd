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
functionality will soon be added for Markov  (i.e., "clock-forward") models.

}
\section{Fields}{

\describe{
\item{\code{trans_model}}{Model for health state transitions. Must be an object of class \code{\link{IndivCtstmTrans}}.}

\item{\code{utility_model}}{The model used to predict utility by health state. Must be an object of class \code{\link{StateVals}}.}

\item{\code{cost_models}}{The models used to predict costs by health state. Must be a list of objects of class \code{\link{StateVals}}, where each element of the list represents a different cost category.}

\item{\code{disprog_}}{An object of class "indiv_ctstm_disprog" simulated using \code{sim_disease}. See \link{IndivCtstmTrans} for more details.}

\item{\code{stateprobs_}}{A tidy \code{\link[data.table]{data.table}} of health state probabilities as a function of time simulated using \code{sim_stateprobs}. The the description of the output of \code{sim_stateprobs} in \code{\link{IndivCtstmTrans}}.

}

\item{\code{qalys_}}{A tidy \code{\link[data.table]{data.table}} of quality-adjusted life-years (QALYs) simulated using \code{sim_qalys}. Columns are:

\describe{
\item{sample}{A random sample from the probabilistic sensitivity analysis (PSA).}
\item{strategy_id}{The treatment strategy ID.}
\item{patient_id}{A patient ID.}
\item{dr}{The discount rate.}
\item{qalys}{Simulated QALYs.}
}

}

\item{\code{costs_}}{A tidy \code{\link[data.table]{data.table}} of costs by category simulated using \code{sim_costs}. Columns are:

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

\item{\code{new(trans_model = NULL, disprog = NULL, utility_model = NULL, cost_models = NULL)}}{

  Constructor for the class. \code{trans_model} must be specified in order to simulate disease progression with \code{sim_disease}. However, if an object of class "indiv_ctstm_disprog" has already been simulated using \code{\link{IndivCtstmTrans}}, then \code{disprog_} can be initialized using \code{disprog}.  
  
}
  
  \item{\code{sim_disease(max_t = 100, max_age = 100)}}{
    Simulate disease progression.
    
   \itemize{
     \item \code{max_t}: Maximum number of time periods to simulate the model.
    \item \code{max_age}: Maximum age to simulate each patient until. Default is 100.
   }
     
    The simulated output is stored in \code{disprog_}. 
  }
  
  \item{\code{sim_stateprobs(t)}}{
    Simulate the probability of being in each health state as a function of time using the simulation output stored in \code{disprog_}.
    
     \itemize{
     \item \code{t}: A numeric vector of times.
     }
     
    The simulated output is stored in \code{stateprobs_}.  
  }
  
  \item{\code{sim_qalys(dr = .03, type = c("predict", "random"), by_patient = TRUE, stat = mean)}}{
    Compute simulated (discounted) QALYs using the simulation output stored in \code{disprog_} by random sample from the PSA, treatment strategy, health state, and 
    (optionally) patient.
    
     \itemize{
     \item \code{dr}: Discount rate to apply to QALYs. May be a vector in which case QALYs are calculated for each element in \code{dr}.
     \item \code{type}: \code{predict} for mean values or \code{random} for random samples as in \code{$sim()} in \code{\link{StateVals}}.
     \item \code{by_patient}: If \code{TRUE}, QALYs are computed at the patient level. If \code{FALSE}, QALYs are averaged across patients by
     health state.
    \item \code{stat}: When \code{by_patient = FALSE}, the statistic used to summarize QALYs by health state. The default is the mean but any R function can be passed to \code{stat}.    
     }
     
     The simulated output is stored in \code{qalys_}.
  } 
  
  \item{\code{sim_costs(dr = .03, type = c("predict", "random"), by_patient = TRUE, stat = mean, max_t = Inf)}}{
    Compute simulated (discounted) costs using the simulation output stored in \code{disprog_} by random sample from the PSA, treatment strategy, health state, and 
    (optionally) patient.
    
     \itemize{
     \item \code{dr}: Discount rate to apply to costs. May be a vector in which case costs are calculated for each element in \code{dr}.
      \item \code{type}: \code{predict} for mean values or \code{random} for random samples as in \code{$sim()} in \code{\link{StateVals}}.
     \item \code{by_patient}: If \code{TRUE}, costs are computed at the patient level. If \code{FALSE}, costs are averaged across patients by
     health state.      
      \item \code{stat}: When \code{by_patient = FALSE}, the statistic used to summarize costs by health state. The default is the mean but any R function can be passed to \code{stat}.  
      \item \code{max_t}: Maximum time duration to compute costs once a patient has entered a (new) health state. By default, equal to \code{Inf}, 
      so that costs are computed over the entire duration that a patient is in a given health state. If time varies by each cost category, then time can also be passed
      as a numeric vector of length equal to the number of cost categories (e.g., \code{c(1, 2, Inf, 3)} for a model with 4 cost categories).
     }
     
    The simulated output is stored in \code{costs_}. 
  },
  
  \item{\code{summarize(stat = mean)}}{
    Produce a summary measure of costs and QALYs by random sample from the PSA and treatment strategy.
    
     \itemize{
     \item \code{stat}: The statistic used to summarize costs and QALYs. The default is the mean but any R function can be passed to \code{stat}.
     }
     
     Returns an object of summary cost-effectiveness output of class \code{\link{ce}}. 
  }    
  

}

}

\seealso{
\code{\link{create_IndivCtstmTrans}}, \code{\link{IndivCtstmTrans}}
}


\examples{
library("flexsurv")

# Treatment strategies,target population, and health states
dt_strategies <- data.frame(strategy_id = c(1, 2, 3))
dt_patients <- data.frame(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
dt_states <- data.frame(state_id = c(1, 2))

# Parameter estimation
n_samples <- 2

## Multi-state model
fits <- vector(length = 3, mode = "list")
surv_dat <- data.frame(ctstm3_exdata$transitions)
for (i in 1:length(fits)){
  fits[[i]] <- flexsurvreg(Surv(years, status) ~ strategy_id, 
                           data = surv_dat,
                           subset = (trans == i),
                           dist = "weibull")
}
fits <- flexsurvreg_list(fits)
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))
              
## Utility
beta_params <- mom_beta(ctstm3_exdata$utility$mean, ctstm3_exdata$utility$se)
utility1 <- rbeta(n_samples, shape1 = beta_params$shape1[1], 
                  shape2 = beta_params$shape2[1])
utility2 <- rbeta(n_samples, shape1 = beta_params$shape1[2], 
                  shape2 = beta_params$shape2[2])
utility_params <- params_lm(coefs  = cbind(utility1, utility2))

# Costs
drugcost_params <- params_lm(matrix(rep(ctstm3_exdata$costs$drugs$costs, n_samples), 
                             nrow = n_samples, byrow = TRUE))

# Decision model
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients,
                        states = dt_states)
                        
## Create transition model (i.e., multi-state model)
transmod_data <- expand_hesim_data(hesim_dat)
transmod <- create_IndivCtstmTrans(fits, data = transmod_data, trans_mat = tmat,
                                   n = n_samples)
                        
## Create costs and utility models
statevals_edata <- expand_hesim_data(hesim_dat, by = c("strategies", "patients",
                                                       "states"))
input_dat <- create_input_data(formula_list(mu = ~ -1 + factor(state_id)), 
                               data = statevals_edata)
utilmod <- StateVals$new(data = input_dat, params = utility_params)
drugcostmod <- StateVals$new(data = input_dat, params = drugcost_params)

## Simulate outcomes
ictstm <- IndivCtstm$new(trans_model = transmod,
                         utility_model = utilmod,
                         cost_models = list(drugs = drugcostmod))
head(ictstm$sim_disease()$disprog_$sim)
head(ictstm$sim_stateprobs(t = c(0, 5, 10))$stateprobs_[t == 5])
ictstm$sim_qalys(dr = .03)
ictstm$sim_costs(dr = .03)
head(ictstm$summarize())
}