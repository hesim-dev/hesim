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

#' Simulate state probabilities from survival curves
#'
#' Simulate health state probabilities from a [`survival`] object using partitioned
#' survival analysis.
#' 
#' @param x An object of class [`survival`].
#' @param ... Further arguments passed to or from other methods.
#' @return 
#' A `stateprobs` object.
#' @examples 
#' library("data.table")
#' library("survival")
#' 
#' # This example shows how to simulate a partitioned survival model by
#' # manually constructing a "survival" object. We will consider a case in which
#' # Cox proportional hazards models from the survival package---which are not
#' # integrated with hesim---are used for parameter estimation. We will use 
#' # point estimates in the example, but bootstrapping, Bayesian modeling,
#' # or other techniques could be used to draw samples for a probabilistic 
#' # sensitivity analysis. 
#' 
#' # (0) We first setup our model per usual by defining the treatment strategies,
#' # target population, and health states
#' hesim_dat <- hesim_data(
#'   strategies = data.table(strategy_id = 1:3,
#'                           strategy_name = c("SOC", "New 1", "New 2")),
#'   patients = data.table(patient_id = 1:2,
#'                         female = c(0, 1),
#'                         grp_id = 1),
#'   states = data.table(state_id = 1:2,
#'                       state_name = c("Stable", "Progression"))
#' )
#'
#' # (1) Next we will estimate Cox models with survival::coxph(). We illustrate 
#' # by predicting progression free survival (PFS) and overall survival (OS)
#' ## Fit models
#' onc3_pfs_os <- as_pfs_os(onc3, patient_vars = c("patient_id", "female",
#'                                                 "strategy_name"))
#' fit_pfs <- coxph(Surv(pfs_time, pfs_status) ~ strategy_name + female,
#'                  data = onc3_pfs_os)
#' fit_os <- coxph(Surv(os_time, pfs_status) ~ strategy_name + female,
#'                 data = onc3_pfs_os)
#'
#' ## Predict survival on input data
#' surv_input_data <- expand(hesim_dat)
#' times <- seq(0, 14, 1/12)
#' predict_survival <- function(object, newdata, times) {
#'   surv <- summary(survfit(object, newdata = newdata, se.fit = FALSE),
#'                   t = times)
#'   pred <- newdata[rep(seq_len(nrow(newdata)), each = length(times)), ]
#'   pred[, sample := 1] # Point estimates only in this example
#'   pred[, time := rep(surv$time, times = nrow(newdata))]
#'   pred[, survival := c(surv$surv)]
#'   return(pred[, ])
#' }
#' pfs <- predict_survival(fit_pfs, newdata = surv_input_data, times = times)
#' os <- predict_survival(fit_os, newdata = surv_input_data, times = times)
#' surv <- rbind(
#'   as.data.table(pfs)[, curve := 1L],
#'   as.data.table(os)[, curve := 2L]
#' )
#'
#' ## Convert predictions to a survival object
#' surv <- survival(surv, t = "time")
#' \dontrun{autoplot(surv)}
#' 
#' # (2) We can then compute state probabilities from the survival object
#' stprobs <- sim_stateprobs(surv)
#' 
#' # (3) Finally, we can use the state probabilities to compute QALYs and costs
#' ## A dummy utility model to illustrate
#' utility_tbl <- stateval_tbl(
#'   data.table(state_id = 1:2,
#'              est = c(1, 1)
#'   ),
#'   dist = "fixed"
#' )
#' utilitymod <- create_StateVals(utility_tbl, 
#'                                hesim_data = hesim_dat,
#'                                n = 1)
#'
#' ## Instantiate Psm class and compute QALYs
#' psm <- Psm$new(utility_model = utilitymod)
#' psm$stateprobs_ <- stprobs
#' psm$sim_qalys()
#' psm$qalys_
#' 
#' @seealso [`survival`]
#' @export
sim_stateprobs.survival <- function(x, ...) {
  state_id <- NULL
  
  # Size
  if (is.null(attr(x, "size"))) stop("'size' attribute missing from 'x'.")
  n_samples <- attr(x, "size")[["n_samples"]]
  n_strategies <- attr(x, "size")[["n_strategies"]]
  n_patients <- attr(x, "size")[["n_patients"]]
  n_states <- attr(x, "size")[["n_states"]]
  n_times <- attr(x, "size")[["n_times"]]
  
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
  if (!"patient_wt" %in% colnames(x)) stprobs[, ("patient_wt") := NULL]
  setattr(stprobs, "class", 
          c("stateprobs", "data.table", "data.frame"))
  setattr(stprobs, "size", 
          c(n_samples = n_samples, n_strategies = n_strategies, 
            n_patients = n_patients, n_states = n_states, n_times = n_times))
  return(stprobs[, ])
}