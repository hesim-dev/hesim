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
\item{\code{params}}{An object of class \code{\link{tparams_transprobs}}.}

\item{\code{input_mats}}{An object of class \code{\link{input_mats}}.}

\item{\code{start_stateprobs}}{A vector with length equal to the number of health states containing the probability that the cohort is in each health state at the start of the simulation. For example, if there were three states and the cohort began the simulation in state 1 then \code{start_stateprobs = c(1, 0, 0)}. If \code{NULL}, then a vector with the first element equal to 1 and all remaining elements equal to 0. }
}
}

\section{Methods}{

\describe{

\item{\code{new(params, input_mats = NULL, start_stateprobs = NULL)}}{

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
# Example to come
}