#' Coerce to `data.table`
#' 
#' Creates a `data.table` that combines the transition probability matrices 
#' and ID columns from a [tparams_transprobs] object. This is often useful for 
#' debugging. 
#' @param x A [tparams_transprobs] object.
#' 
#' @seealso [tparams_transprobs]
#' @return A `data.table` with one row for each transition probability matrix.
#' @export
as.data.table.tparams_transprobs <- function(x){
  probs <- as_tpmatrix(x$value)
  colnames(probs) <- paste0("prob_", 1:ncol(probs)) 
  id_dt <- as.data.table(x[c("sample", "strategy_id", "patient_id")])
  time_dt <- x$time_intervals[match(x$time_id, x$time_intervals$time_id)]
  x_dt <- data.table(id_dt, time_dt, probs)
  for (v in c("n_samples", "n_strategies", "n_patients", "n_times")){
    setattr(x_dt, v, x[[v]])
  }
  return(x_dt)
}