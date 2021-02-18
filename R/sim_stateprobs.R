#' Simulated state probabilities
#'
#' A generic function to simulate state probabilities and create an object of 
#' class [`stateprobs`].
#' 
#' @param x An object of the appropriate class.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' A [`stateprobs`] object.
#' @export
#' @seealso [`sim_stateprobs.survival`]
sim_stateprobs <- function(x, ...) {
  UseMethod("sim_stateprobs")
}

#' State probabilities from survival curves
#'
#' Simulate health state probabilities from a [`survival`] object using partitioned
#' survival analysis.
#' 
#' @param x An object of class [`survival`].
#' @param ... Further arguments passed to or from other methods.
#' @return 
#' A `stateprobs` object.
#' @export
sim_stateprobs.survival <- function(x, ...) {
  # Size
  n_samples <- length(unique(x$sample))
  n_strategies <- length(unique(x$strategy_id))
  n_patients <- length(unique(x$patient_id))
  n_states <- length(unique(x$curve)) + 1L
  n_times <- length(unique(x$t))
  
  # Simulate state probabilities
  res <- C_psm_sim_stateprobs(x,
                              n_samples = n_samples,
                              n_strategies = n_strategies,
                              n_patients = n_patients,
                              n_states = n_states,
                              n_times = n_times)
  prop_cross <- res$n_crossings/nrow(res$stateprobs)
  if (prop_cross > 0){
    warning(paste0("Survival curves crossed ", round(prop_cross * 100, 2), 
                   " percent of the time."),
            call. = FALSE)
  }
  
  # Create object and set attributes
  stprobs <- data.table(res$stateprobs)
  stprobs[, state_id := state_id + 1]
  stprobs[, sample := sample + 1]
  check_patient_wt(x, stprobs)
  setattr(stprobs, "class", 
          c("stateprobs", "data.table", "data.frame"))
  setattr(stprobs, "size", 
          c(n_samples = n_samples, n_strategies = n_strategies, 
            n_patients = n_patients, n_states = n_states, n_times = n_times))
  return(stprobs[, ])
}
