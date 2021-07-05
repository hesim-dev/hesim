#' @importFrom data.table as.data.table
#' @export
data.table::as.data.table

#' Coerce to `data.table`
#' 
#' Creates a `data.table` that combines the transition probability matrices 
#' and ID columns from a [tparams_transprobs] object. This is often useful for 
#' debugging. 
#' @param x A [tparams_transprobs] object.
#' @param prefix,sep Arguments passed to [tpmatrix_names()] for naming
#' the transition probability columns. The `states` argument is based on
#' the column names (i.e., names of the second dimension) of the `$value`
#' element of `x`; if `NULL`, then states are named `s1`, ..., `sh` where h is 
#' the number of states.
#' 
#' @seealso [tparams_transprobs]
#' @return A `data.table` with one row for each transition probability matrix.
#' 
#' @examples 
#' # Create tparams_transprobs object
#' hesim_dat <- hesim_data(strategies = data.frame(strategy_id = 1:2),
#'                         patients = data.frame(patient_id = 1:3))
#' input_data <- expand(hesim_dat, by = c("strategies", "patients"))    
#' tpmat_id <- tpmatrix_id(input_data, n_samples = 2)      
#' p_12 <- runif(nrow(tpmat_id), .6, .7) + 
#'   .05 * (tpmat_id$strategy_id == 2)
#' tpmat <- tpmatrix(
#'   C, p_12,
#'   0, 1
#' )
#' tprobs <- tparams_transprobs(tpmat, tpmat_id)
#' 
#' # Convert to data.table
#' as.data.table(tprobs)
#' as.data.table(tprobs, prefix = "")
#' as.data.table(tprobs, prefix = "", sep = ".")
#' 
#' @export
as.data.table.tparams_transprobs <- function(x, prefix = "prob_", sep = "_"){
  probs <- as_tbl2(x$value, prefix = prefix, sep = sep)
  id_dt <- as.data.table(x[c("sample", "strategy_id", "patient_id")])
  time_dt <- x$time_intervals[match(x$time_id, x$time_intervals$time_id)]
  x_dt <- data.table(id_dt, time_dt, probs)
  for (v in c("n_samples", "n_strategies", "n_patients", "n_times")){
    setattr(x_dt, v, x[[v]])
  }
  return(x_dt)
}