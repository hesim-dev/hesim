% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tparams_transprobs.R
\name{tparams_transprobs}
\alias{tparams_transprobs}
\alias{print.tparams_transprobs}
\alias{tparams_transprobs.array}
\alias{tparams_transprobs.data.table}
\alias{tparams_transprobs.data.frame}
\alias{tparams_transprobs.tpmatrix}
\title{Transition probabilities}
\usage{
tparams_transprobs(object, ...)

\method{tparams_transprobs}{array}(
  object,
  tpmatrix_id = NULL,
  times = NULL,
  grp_id = NULL,
  patient_wt = NULL,
  ...
)

\method{tparams_transprobs}{data.table}(object, ...)

\method{tparams_transprobs}{data.frame}(object, ...)

\method{tparams_transprobs}{tpmatrix}(object, tpmatrix_id, ...)
}
\arguments{
\item{object}{An object of the appropriate class.}

\item{...}{Further arguments passed to or from other methods. Currently unused.}

\item{tpmatrix_id}{An object of class \code{\link{tpmatrix_id}} (or an equivalent
\code{data.table} with the same ID columns as returned by \code{tpmatrix_id()}).}

\item{times}{An optional numeric vector of distinct times to pass to
\link{time_intervals} representing time intervals indexed by the 4th dimension of
the array. May either be the start or the end of intervals.
This argument is not required if there is only one time interval.}

\item{grp_id}{An optional numeric vector of integers denoting the subgroups. Must
be the same length as the 3rd dimension of the array.}

\item{patient_wt}{An optional numeric vector denoting the weight to apply to each
patient within a subgroup. Must be the same length as the 3rd dimension of the array.}
}
\value{
An object of class \code{tparams_transprobs},
which is a list containing \code{value} and relevant ID attributes. The element
\code{value} is an array of predicted transition probability matrices from the
probability distribution of the underlying statistical model. Each matrix in
\code{value} is a prediction for a \code{sample}, \code{strategy_id}, \code{patient_id}, and
optionally \code{time_id} combination.
}
\description{
Create a list containing predicted transition probabilities at discrete times.
Since the transition probabilities have presumably already been predicted
based on covariate values, no input data is required for
simulation. The class can be instantiated from either an \code{array},
a \code{data.table}, a \code{data.frame}, or a \code{\link{tpmatrix}}. This is the object in
\code{hesim} used to specify the transition probabilities required to simulate
Markov chains with the \code{\link{CohortDtstmTrans}} class.
}
\details{
The format of \code{object} depends on its class:
\describe{
\item{array}{
Either a 3D or a 6D array is possible.
\itemize{
\item If a 3D array, then each slice is a
square transition probability matrix. In this case
\code{tpmatrix_id} is required and each matrix slice corresponds to the same
numbered row in \code{tpmatrix_id}. The number of matrix slices must equal the number
of rows in \code{tpmatrix_id}.

\item If a 6D array, then the dimensions of the array should be indexed as follows:
1st (\code{sample}), 2nd (\code{strategy_id}), 3rd (\code{patient_id}),
4th (\code{time_id}), 5th (rows of transition matrix), and
6th (columns of transition matrix). In other words, an index of
\verb{[s, k, i, t]} represents the transition matrix for the \code{s}th
sample, \code{k}th treatment strategy, \code{i}th patient, and \code{t}th
time interval.
}
}

\item{data.table}{Must contain the following:
\itemize{
\item ID columns for the parameter sample (\code{sample}),
treatment strategy (\code{strategy_id}), and patient (\code{patient_id}).
If the number of time intervals is greater than 1 it must also contain the
column \code{time_start} denoting the starting time of a time interval. A column
\code{patient_wt} may also be used to denote the weight to apply to each
patient.
\item Columns for each element of the transition probability matrix.
They should be prefixed with "prob_" and ordered rowwise.
For example, the following columns would be used for a 2x2 transition
probability matrix:
\code{prob_1} (1st row, 1st column),
\code{prob_2} (1st row, 2nd column),
\code{prob_3} (2nd row, 1st column), and
\code{prob_4} (2nd row, 2nd column).
}
}

\item{data.frame}{Same as \code{data.table}.}

\item{tpmatrix}{An object of class \code{\link{tpmatrix}}.}
}

A \code{tparams_transprobs} object is also instantiated when creating a
cohort discrete time state transition model using \code{\link[=define_model]{define_model()}}.
}
\examples{
hesim_dat <- hesim_data(strategies = data.frame(strategy_id = 1:2),
                        patients = data.frame(patient_id = 1:3))
input_data <- expand(hesim_dat, by = c("strategies", "patients"))    

# tpmatrix objects provide a convenient way to construct
# tparams_transprobs() objects
tpmat_id <- tpmatrix_id(input_data, n_samples = 2)      
p_12 <- runif(nrow(tpmat_id), .6, .7) + 
  .05 * (tpmat_id$strategy_id == 2)
tpmat <- tpmatrix(
  C, p_12,
  0, 1
)
tprobs <- tparams_transprobs(tpmat, tpmat_id)
names(tprobs) # Names of list elements

# Convert to data.table in wide format
as.data.table(tprobs)

# Convert to data.table in long format
as.data.table(tprobs, long = TRUE)

# Summary where each column is a vector
summary(tprobs)
summary(tprobs, probs = c(.025, .975))

# Summary where each column is a matrix
ps <- summary(tprobs, id = tpmat_id, unflatten = TRUE)
ps
ps$mean
}
\seealso{
A \code{tparams_transprobs} object is used to store the "parameters" of
the transition component of a cohort discrete time state transition
model (cDTSTM). You can create such an object with \code{CohortDtstmTran$new()}.

\code{\link[=tpmatrix]{tpmatrix()}} and \code{\link[=tpmatrix_id]{tpmatrix_id()}} provide a convenient way to construct a
\code{tparams_transprobs} object in a flexible way. \code{\link[=define_model]{define_model()}} is, in turn,
a convenient way to construct a \code{\link{tpmatrix}} object using mathematical
expressions; in this case, an entire cDTSTM can be instantiated from a model
definition using \code{\link[=create_CohortDtstm.model_def]{create_CohortDtstm.model_def()}}. Detailed examples
are provided in \code{vignette("markov-cohort")} and
\code{vignette("markov-inhomogeneous-cohort")}

The output of a \code{tparams_transprobs} object is rather verbose. It can be
helpful to check the output by converting it to a \code{data.table} (containing
both the ID variables and flattened transition probability matrices)
with \code{\link[=as.data.table.tparams_transprobs]{as.data.table.tparams_transprobs()}}. Transition probabilities can
also be summarized (across parameter samples) using
\code{\link[=summary.tparams_transprobs]{summary.tparams_transprobs()}}.
}
