\docType{data}
\name{CohortDtstmTrans}
\alias{CohortDtstmTrans}

\title{Transitions for a cohort discrete time state transition model}
\format{\code{\link{R6Class}} object.}
\usage{
CohortDtstmTrans
}
\description{
Simulate health state transitions in a cohort discrete time state transition model.

}
\section{Fields}{

\describe{
\item{\code{input_mats}}{ An object of class \code{\link{input_mats}}. Since \code{params} currently only supports transition matrices that have already been predicted (i.e., they have been stored in \code{\link{transprob_tbl}}), input matrices are not needed for prediction and \code{input_mats} only needs to contain ID variables indexing the rows of \code{params}.}

\item{\code{params}}{An array of transition matrices where each matrix contains the transition probabilities for a given sample from the probabilistic sensitivity analysis (PSA), treatment strategy, patient, and time interval. Must be sorted by PSA sample, treatment strategy, patient, and time interval.}

\item{\code{start_stateprobs}}{A vector with length equal to the number of health states containing the probability that the cohort is in each health state at the start of the simulation. For example, if there were three states and the cohort began the simulation in state 1 then \code{start_stateprobs = c(1, 0, 0)}.}

}
}

\section{Methods}{

\describe{

\item{\code{new(input_mats, params, start_stateprobs)}}{

  Constructor for the class.
  
}
  
  \item{\code{sim_stateprobs(n_cycles)}}{
    Simulate probability of being in each health state during each model cycle. 
    
     \itemize{
     \item \code{n_cycles}: The number of model cycles to simulate the model for.
     }
  }
  
  }  
  
}

\seealso{
\code{\link{create_CohortDtstmTrans}}
}

\examples{
#' # Observed transitions 
#' # Treatment strategies 1 and 2
#' trans1 <- c(200, 300, 500,
#'             0, 400, 100,
#'             0, 0, 600) 
#' trans2 <- c(300, 300, 400,
#'             0, 550, 50,
#'             0, 0, 450)
#' n_trans <- length(trans1)

#' # Sample transition probabilities by treatment strategy
#' # from a dirichlet distributions
#' transprobs <- transprob_tbl(data.frame(strategy_id = rep(1:2, each = n_trans),
#'                                        patient_id = 1,
#'                                        transition_id = rep(1:n_trans, 2),
#'                                        alpha = c(trans1,trans2)),
#'                             dist = "dirichlet",
#'                             hesim_data = hesim_dat)
#' print(transprobs)
#' create_CohortDtstmTrans(transprobs, n = 3)
#' print(transprobs$params)
}