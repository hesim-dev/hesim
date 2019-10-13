\docType{data}
\name{CohortDtstm}
\alias{CohortDtstm}

\title{Cohort discrete time state transition model}
\format{\code{\link{R6Class}} object.}
\usage{
CohortDtstm
}
\description{
Simulate outcomes from a cohort discrete time state transition model.

}
\section{Fields}{

\describe{
\item{\code{trans_model}}{Model for health state transitions. Must be an object of class \code{\link{CohortDtstmTrans}}.}

\item{\code{utility_model}}{The model used to predict utility by health state. Must be an object of class \code{\link{StateVals}}.}

\item{\code{cost_models}}{The models used to predict costs by health state. Must be a list of objects of class \code{\link{StateVals}}, where each element of the list represents a different cost category.}

\item{\code{stateprobs_}}{An object of class \code{\link{stateprobs}}.}

\item{\code{costs_}}{An object of class \code{costs}.}

\item{\code{qalys_}}{An object of class \code{qalys}.}

}
}

\section{Methods}{

\describe{

\item{\code{new(trans_model = NULL, utility_model = NULL, cost_models = NULL)}}{

  Constructor for the class.
  
}

\item{\code{sim_qalys(dr = .03, method = c("trapz", "riemann_left", "riemann_right"),
                      lys = TRUE}}{
Simulate QALYS using \code{\link{sim_qalys}}.

} 

\item{\code{sim_costs(dr = .03, method = c("trapz", "riemann_left", "riemann_right")}}{

Simulate costs using \code{\link{sim_costs}}.

}  
  
}
}

\seealso{
\code{\link{create_CohortDtstmTrans}}, \code{\link{CohortDtstmTrans}}
}

\examples{
# Example to come
}