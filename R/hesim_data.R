# hesim data -------------------------------------------------------------------

#' Create a data table of treatment lines
#' 
#' Convert a list of treatment lines for multiple treatment strategies to a 
#' `data.table`.
#' @param strategy_list A list where each element is a treatment strategy 
#' consisting of a vector of treatments. 
#' @param strategy_ids A numeric vector denoting the numeric id of each strategy
#' in \code{strategy_list}.
#' @return Returns a `data.table` in tidy format with three columns:
#' \describe{
#' \item{strategy_id}{Treatment strategy ids.}
#' \item{line}{Line of therapy.}
#' \item{treatment_id}{Treatment ID for treatment used at a given line of therapy 
#' within a treatment strategy.}
#' }
#' @examples 
#' strategies <- list(c(1, 2, 3),
#'                   c(1, 2))
#' create_lines_dt(strategies)
#' @export
create_lines_dt <- function(strategy_list, strategy_ids = NULL){
  treatments <- unlist(strategy_list, use.names = FALSE)
  if(!is.numeric(treatments)){
    stop("Elements in 'strategy_list' should be integers.")
  }
  n_treatments <- unlist(lapply(strategy_list, length), 
                         use.names = FALSE)
  if (is.null(strategy_ids)){
    strategy_ids <- seq_len(length(strategy_list))
  }
  strategies <- rep(strategy_ids, times = n_treatments)
  lines <- unlist(lapply(n_treatments, seq_len), use.names = FALSE)
  return(data.table(strategy_id = strategies,
                    line = lines,
                    treatment_id = treatments))
}

#' Create a data table of health state transitions
#' 
#' Create a data table of health state transitions from a transition matrix describing 
#' the states and transitions in a multi-state model suitable for use with [`hesim_data`]. 
#' @param trans_mat A transition matrix in the format from the [`mstate`][mstate::mstate] package. 
#' See [`IndivCtstmTrans`].
#' @return Returns a [`data.table`] in tidy format with three columns:
#' \describe{
#' \item{transition_id}{Health state transition ID.}
#' \item{from}{The starting health state.}
#' \item{to}{The health state that will be transitions to.}
#' }
#' @examples 
#' tmat <- rbind(c(NA, 1, 2),
#'               c(NA, NA, 3),
#'               c(NA, NA, NA))
#' create_trans_dt(tmat)
#' @export
create_trans_dt <- function(trans_mat){
  n_rows <- nrow(trans_mat)
  id <- to <- from <- vector(mode = "list", n_rows)
  for (i in 1:n_rows){
    id_i <- trans_mat[i, ]
    id[[i]] <- id_i[!is.na(id_i)]
    from[[i]] <- rep(i, length(id[[i]]))
    names(from[[i]]) <- rep(rownames(trans_mat)[i], length(from[[i]]))
    to[[i]] <- which(!is.na(id_i))
  }
  id <- do.call("c", id)
  from <- do.call("c", from)
  to <- do.call("c", to)
  x <- data.table(transition_id = id,
                  from = from,
                  to = to)
  if (!is.null(names(from)) & !is.null(names(to))){
    x$from_name <- names(from)
    x$to_name <- names(to)
  }
  setorderv(x, "transition_id")
  return(x)
}

#' Data for health economic simulation modeling
#' 
#' A list of tables required for health economic simulation modeling. This object
#' is used to setup models by defining the treatment strategies, target population,
#' and model structure.
#' @param strategies A table of treatment strategies. Must contain the column 
#' `strategy_id` denoting a unique strategy. Other columns are variables
#' describing the characteristics of a treatment strategy. 
#' @param patients A table of patients. Must contain the column `patient_id` denoting 
#' a unique patient. The number of rows should be equal to the number of patients 
#' in the model. The table may also include columns for `grp_id` for subgroups and 
#' `patient_wt` specifying the weight to apply to each patient (within a subgroup). 
#' If `grp_id` is `NULL`, then it is assumed that there is only one subgroup. If
#' `patient_wt` is `NULL`. then each patient is given the same weight. Weights 
#' cannot be used in individual-level models because each patient should be
#' weighted equally; that is, weights can only be specified in cohort models.
#' Weights within subgroups are normalized to sum to one. Other columns are 
#' variables describing the characteristics of a patient.
#' @param states A table of health states. Must contain the column
#' `state_id`, which denotes a unique health state. The number of rows should
#' be equal to the number of health states in the model. Other columns can describe the
#' characteristics of a health state.
#' @param transitions A table of health state transitions. Must contain the column
#' `transition_id`, which denotes a unique transition; `from`, which denotes
#' the starting health state; and `to` which denotes the state that will be
#' transitioned to.
#' 
#' @note Each table must either be a `data.frame` or `data.table`. All ID variables 
#' within each table must be numeric vectors of integers and should be of the form
#' 1,2,...N where N is the number of unique values of the ID variable. 
#' 
#' @return Returns an object of class `hesim_data`, which is a list of data tables for
#' health economic simulation modeling.
#' @seealso [`expand.hesim_data()`], [`get_labels()`]
#' @examples 
#' strategies <- data.frame(strategy_id = c(1, 2))
#' patients <- data.frame(patient_id = seq(1, 3), age = c(65, 50, 75),
#'                        gender = c("Female", "Female", "Male"))
#' states <- data.frame(state_id =  seq(1, 3),
#'                      state_var = c(2, 1, 9))
#' hesim_dat <- hesim_data(strategies = strategies,
#'                         patients = patients,
#'                         states = states)
#' @export
hesim_data <- function(strategies, patients, states = NULL,
                       transitions = NULL){
  object <- new_hesim_data(strategies = strategies, patients = patients,
                           states = states, transitions = transitions)
  return(check(object))
}

new_hesim_data <- function(strategies, patients, states = NULL,
                           transitions = NULL){
  data <- list()
  data$strategies <- strategies 
  data$patients <- patients
  data$states <- states
  data$transitions <- transitions
  class(data) <- c("hesim_data")
  return(data)
}

check.hesim_data <- function(x){
  # Strategies
  if (is.null(x$strategies)){
    stop("'strategies' cannot be NULL'.",
         call. = FALSE)
  }
  check_hesim_data_type(x$strategies, "strategies")
  if (!"strategy_id" %in% colnames(x$strategies)){
    stop("'strategies' must contain the column 'strategy_id'.",
         call. = FALSE)
  }
  
  # Patients
  if (!is.null(x$patients)){
    check_hesim_data_type(x$patients, "patients")
    if (!"patient_id" %in% colnames(x$patients)){
      stop("'patients' must contain the column 'patient_id'.", 
           call. = FALSE)
    }
  }
  
  # States
  if (!is.null(x$states)){
    check_hesim_data_type(x$states, "states")
    if (!"state_id" %in% colnames(x$states)){
      stop("'states' must contain the column 'state_id'.", 
           call. = FALSE)
    }
  }
  
  # Transitions
  if (!is.null(x$transitions)){
    check_hesim_data_type(x$transitions, "transitions")
    if (!"transition_id" %in% colnames(x$transitions)){
      stop("'transitions' must contain the column 'transition_id'.", 
           call. = FALSE)
    }
  }
  
  return(x)
}

check_hesim_data_type <- function(tbl, tbl_name){
  if(!is.data.frame(tbl) & !is.data.table(tbl)){
    msg <- paste0("'", tbl_name, "'", " must be a data.frame or data.table.")
    stop(msg, call. = FALSE)
  }
}

#' Expand object
#' 
#' A generic function for "expanding" an object. Only used for 
#' `hesim_data` objects with [expand.hesim_data()].
#' @export
#' @keywords internal
#' @seealso [expand.hesim_data()]
expand <- function(object, by, times){
  UseMethod("expand", object)
}

#' Expand hesim_data
#' 
#' Create a data table in long format from all combinations of specified tables 
#' from an object of class [hesim_data] and optionally time intervals. See "Details" for 
#' an explanation of how the expansion is done.
#' @param object An object of class `hesim_data`.
#' @param by A character vector of the names of the data tables in `hesim_data` to expand by.
#' @param times Either a numeric vector of distinct times denoting the start of time intervals or
#' a [time_intervals] object. 
#' @details This function is similar to [expand.grid()], but works for data frames or data tables. 
#' Specifically, it creates a `data.table` from all combinations of the supplied tables in `object`
#' and optionally the start of times intervals in `times`. 
#' The supplied tables are determined using the `by` argument. The resulting dataset is sorted by 
#' prioritizing ID variables as follows: (i) `strategy_id`, (ii) `patient_id`,
#' (iii) the health-related ID variable (either `state_id` or `transition_id`, and
#' (iv) the time intervals from `times`.
#' @return An object of class `expanded_hesim_data`, which is a `data.table` with an "id_vars" 
#' attribute containing the names of the ID variables in the data table and, if `times` is 
#' not `NULL`, a `time_intervals` object derived from `times`.
#' @examples 
#' strategies <- data.frame(strategy_id = c(1, 2))
#' patients <- data.frame(patient_id = seq(1, 3), age = c(65, 50, 75),
#'                           gender = c("Female", "Female", "Male"))
#' states <- data.frame(state_id =  seq(1, 3),
#'                      state_var = c(2, 1, 9))
#' hesim_dat <- hesim_data(strategies = strategies,
#'                         patients = patients,
#'                         states = states)
#' expand(hesim_dat, by = c("strategies", "patients"))
#' expand(hesim_dat, by = c("strategies", "patients"),
#'        times = c(0, 2, 10))
#' @export
expand.hesim_data <- function(object, by = c("strategies", "patients"),
                              times = NULL){
  if ("transitions" %in% by & "states" %in% by){
    stop("Cannot expand by both 'transitions' and 'states'.", call. = FALSE)
  }
  allowed_tables <- names(hesim_data_sorting_map()[
      !names(hesim_data_sorting_map()) == "times"
    ])
  if (!all(by %in% allowed_tables)){
    stop("One of the elements specified in 'by' is not a table in 'hesim_data'.",
         call. = FALSE)
  }
  sorted_by <- hesim_data_sorted_by(by)
  tbl_list <- object[sorted_by]
  if (!is.null(times)){
    if (inherits(times, "time_intervals")){
      tintervals <- times
    } else{
      tintervals <- time_intervals(times)
    }
    sorted_by <- hesim_data_sorted_by(c(by, "times"))
    tbl_list <- c(tbl_list, list(times = tintervals))
  }
  for (i in 1:length(tbl_list)){
    if (is.null(tbl_list[[i]])){
      stop("Cannot merge a NULL data table.")
    }
  }
  tbl_list <- lapply(tbl_list, function(x) as.data.frame(x))
  dat <- data.table(Reduce(function(...) merge(..., by = NULL, sort = FALSE),
                           tbl_list))
  sort_hesim_data(dat, sorted_by)
  id_cols <- unlist(hesim_data_sorting_map()[sorted_by])
  nonid_cols <- colnames(dat)[!colnames(dat) %in% id_cols]
  dat <- dat[, c(id_cols, nonid_cols), with = FALSE]
  setattr(dat, "id_vars", unname(id_cols))
  if (!is.null(times))  setattr(dat, "time_intervals", tintervals)
  setattr(dat, "class", c("expanded_hesim_data", "data.table", "data.frame"))
  return(dat)
}

#' @export
`[.expanded_hesim_data` <- function(x, ...) {
  out <- NextMethod()
  setattr(out, "id_vars", attributes(x)$id_vars)
  return(out)
}

#' @export
cbind.expanded_hesim_data <- function(x, ...) {
  l <- list(x, ...)
  out <- cbind.data.frame(l)
  setattr(out, "class", c("expanded_hesim_data", "data.table", "data.frame"))
  setattr(out, "id_vars", attributes(x)$id_vars)
  return(out)
}

hesim_data_sorting_map <- function(){
  list(strategies = "strategy_id",
       patients = "patient_id",
       transitions = "transition_id",
       states = "state_id",
       times = "time_id")
}

hesim_data_sorted_by <- function(by){
  sorted_map <- hesim_data_sorting_map()
  sorted_by <- names(sorted_map)[names(sorted_map) %in% by]
  return(sorted_by)
}

sort_hesim_data <- function(data, sorted_by){
  setorderv(data, unlist(hesim_data_sorting_map()[sorted_by]))
}

# ID attributes ----------------------------------------------------------------
#' Time intervals
#' 
#' Create a table of time intervals given a vector or data frame of unique times.
#' This would typically be passed to [id_attributes].
#' 
#' @param times Either a vector of starting times for each interval or a
#'  `data.frame` with at least one column named `time_start`.
#' @return An object of class `time_intervals` that inherits from
#'  `data.table` in the same format as `time_intervals` as 
#' described in [id_attributes].
#' @seealso [id_attributes]
#' @examples
#' time_intervals(c(0, 3, 5))
#' time_intervals(data.frame(time_start = c(0, 3, 5),
#'                           time_cat = c("Time <= 3", "3 < Time <= 5", 
#'                                        "Time > 5")))
#' @export
time_intervals <- function(times){
  if (inherits(times, "data.frame")){
      if (!"time_start" %in% colnames(times)){
        stop(paste0("If 'time_start' is a data frame, then 'times' ",
                    "must contain a column named 'time_start'."))
      }
      times <- data.table(times)
      if (!any(times[["time_start"]] <= 0)){
        stop(paste0("If 'times' is a data.frame, then the column ", 
                    "'time_start' must contain at least one value <= 0"))
      }
      if (any(is.infinite(times[["time_start"]]))){
        stop(paste0("If 'times' is a data.frame, then the column ", 
                    "'time_start' cannot contain a value equal to 'Inf'."))
      }
      times[, ("time_start") := as.numeric(get("time_start"))]
      setorderv(times, "time_start")
      time_intervals <- data.table(time_id = 1:nrow(times), 
                                   times)
    } else{
      times <- times[!is.infinite(times)]
      times <- times[times >= 0]
      if (!any(times <= 0)) times <- c(0, times)
      time_intervals <- data.table(time_id = 1:length(times), 
                                   time_start = sort(as.numeric(times)))
    }
  time_intervals[, ("time_stop") := shift(get("time_start"), type = "lead")]
  time_intervals[is.na(get("time_stop")), "time_stop" := Inf]
  setattr(time_intervals, "class", c("time_intervals", "data.table", "data.frame"))
  return(time_intervals[, ])
}

#' Attributes for ID variables
#' 
#' Stores metadata related to the ID variables used to index [input_mats] 
#' and [transformed parameter objects][tparams] already predicted from covariates.
#' 
#' @param strategy_id A numeric vector of integers denoting the treatment strategy.
#' @param n_strategies A scalar denoting the number of unique treatment strategies.
#' @param patient_id A numeric vector of integers denoting the patient.
#' @param n_patients A scalar denoting the number of unique patients.
#' @param state_id A numeric vector of integers denoting the health state.
#' @param n_states A scalar denoting the number of unique health states.
#' @param transition_id A numeric vector denoting the 
#' health state transition. This is only used for state transition models. 
#' @param n_transitions A scalar denoting the number of unique transitions. 
#' @param time_id A numeric vector of integers denoting a unique time interval.
#' @param time_intervals A `data.table` denoting unique time intervals. Must 
#' contain the columns `time_id`, `time_start`, and `time_stop`.
#' `time_start` is the starting time of an interval and `time_stop` is
#' the stopping time of an interval. Following the [survival][survival::tmerge] package,
#' time intervals are closed on the right and
#' open on the left (except in the final interval where `time_stop` is equal to 
#' infinity). 
#' @param n_times A scalar denoting the number of time intervals. Equal to the
#' number of rows in `time_intervals`.
#' @param sample A numeric vector of integer denoting the sample from the posterior
#' distribution of the parameters. 
#' @param n_samples A scalar denoting the number of samples.
#' @param grp_id An optional numeric vector of integers denoting the subgroup. 
#' @param patient_wt An optional numeric vector denoting the weight to apply to each patient
#' within a subgroup. 
#'  
#' @details When using the ID variables to index [input_mats], the sorting order 
#' should be the same as specified in [expand.hesim_data()]; that is,
#' observations must be sorted by prioritizing: (i) `strategy_id`, (ii) `patient_id`, 
#' (iii) the health-related ID variable (either `state_id` or `transition_id`), 
#' and (iv) `time_id`. When using ID variables are used to index transformed parameter 
#' objects and `sample` is used for indexing, then observations must be sorted by
#' prioritizing: (i) `sample`, (ii) `strategy_id`, (iii) `patient_id`,
#' (iv) the health-related ID variable, and (v) `time_id`. 
#'   
#' @seealso [hesim_data()],[expand.hesim_data()], [input_mats]
#' @export
id_attributes <- function(strategy_id, n_strategies,
                          patient_id, n_patients,
                          state_id = NULL, n_states = NULL,
                          transition_id = NULL, n_transitions = NULL,
                          time_id = NULL, time_intervals = NULL, n_times = NULL,
                          sample = NULL, n_samples = NULL,
                          grp_id = NULL, patient_wt = NULL){
  object <- new_id_attributes(strategy_id, n_strategies,
                              patient_id, n_patients,
                              state_id, n_states,
                              transition_id, n_transitions,
                              time_id, time_intervals, n_times,
                              sample, n_samples,
                              grp_id, patient_wt)
  check(object)
  return(object)
}

new_id_attributes <- function(strategy_id, n_strategies,
                              patient_id, n_patients,
                              state_id = NULL, n_states = NULL,
                              transition_id = NULL, n_transitions = NULL,
                              time_id = NULL, time_intervals = NULL, n_times = NULL,
                              sample = NULL, n_samples = NULL,
                              grp_id = NULL, patient_wt = NULL){
  stopifnot(is.numeric(strategy_id))
  stopifnot(is.numeric(n_strategies))
  stopifnot(is.numeric(patient_id))
  stopifnot(is.numeric(n_patients))
  stopifnot(is.numeric(state_id) | is.null(state_id))
  stopifnot(is.numeric(n_states) | is.null(n_states))
  stopifnot(is.numeric(transition_id) | is.null(n_transitions))
  stopifnot(is.numeric(n_transitions) | is.null(n_transitions))
  stopifnot(is.numeric(time_id) | is.null(time_id))
  stopifnot(is.data.table(time_intervals) | is.null(time_intervals))
  stopifnot(is.numeric(n_times) | is.null(n_times))
  stopifnot(is.numeric(sample) | is.null(sample))
  stopifnot(is.numeric(n_samples) | is.null(n_samples))
  stopifnot(is.numeric(grp_id) | is.null(grp_id))
  stopifnot(is.numeric(patient_wt) | is.null(patient_wt))
  if (is.null(grp_id)) grp_id <- rep(1, length(patient_id))
  
  object <- list(strategy_id = strategy_id, n_strategies = n_strategies,
                 patient_id = patient_id, n_patients = n_patients,
                 state_id = state_id, n_states = n_states,
                 transition_id = transition_id, n_transitions = n_transitions,
                 time_id = time_id, time_intervals = time_intervals, n_times = n_times,
                 sample = sample, n_samples = n_samples,
                 grp_id = grp_id, patient_wt = patient_wt)
  object[sapply(object, is.null)] <- NULL
  class(object) <- "id_attributes"
  return(object)
}

#' @rdname check
check.id_attributes <- function(object){
  # ID variables to check
  id_vars <- c("sample", "strategy_id", "patient_id", "state_id", 
               "transition_id", "time_id")
  id_vars_n <- c("n_samples", "n_strategies", "n_patients", "n_states",
                 "n_transitions", "n_times")
  keep <- which(id_vars %in% names(object))
  id_vars <- id_vars[keep]
  id_vars_n <- id_vars_n[keep]
  
  # Helper function
  str_list <- function(v) {
    if (length(v) == 1) {
      return(v)
    } else if (length(v) == 2) {
      return(paste0(v[1], " and ", v[2]))
    } else{
      return(paste0(paste(v[1:(length(v) - 1)], collapse = ", "),
            ", and ", v[length(v)]))
    }
  }
  
  # Check that id variables have correct size
  for (i in 1:length(id_vars)){
    if(length(unique(object[[id_vars[i]]])) != object[id_vars_n[i]]){
      msg <- paste0("The number of unique observations in '", id_vars[i], 
                    "' does not equal '", id_vars_n[i], "'.")
      stop(msg, call. = FALSE)
    } 
  }
  
  # Check that each ID vector is same length
  actual_N <- sapply(object[id_vars], length)
  expected_N <- prod(unlist(object[id_vars_n]))
  if(sum(actual_N != expected_N) > 0) {
    stop(paste0("The length of the ID variables is not consistent with the number ",
                "of unique values of each ID variable."), call. = FALSE)
  }
  
  # Check if id variables are sorted properly 
  indices_df <- data.table(do.call("cbind", object[id_vars]))
  sorted_seq <- seq_len(nrow(indices_df))
  indices_df[, "row_num" := sorted_seq]
  setorderv(indices_df, id_vars)
  if(!all(indices_df$row_num == sorted_seq)){
    msg <- paste0("The ID variables are not sorted correctly. The sort priority of the ",
                  "ID variables must be as follows: ", str_list(id_vars), ".")
    stop(msg, call. = FALSE)
  }
  
  # Check if the number of unique observations is correct within groups
  indices_df[, "row_num" := NULL]
  for (i in 2:ncol(indices_df)){
    dt_by <- colnames(indices_df)[1:(i - 1)]
    col <- colnames(indices_df)[i]
    len <- indices_df[, list(len = length(unique(get(col)))), 
                      by = dt_by]$len
    user_n <- object[[id_vars_n[i]]]
    if (!all(unique(len) == user_n)){
      msg <- paste0("The number of unique ", col, " observations within each ",
                     str_list(dt_by), " group must equal ", id_vars_n[i], ".")
      stop(msg, call. = FALSE)
    }
  }
}

# Get the object containing ID attributes
get_id_object <- function(x){
  if (is.null(x$input_data)){
    return(x$params)
  } else{
    return(x$input_data)
  }
}

# Get sizes from an ID object
get_size <- function(x) {
  y <- get_id_object(x)
  return(c(
    n_samples = get_n_samples(x$params),
    n_strategies = y$n_strategies,
    n_patients = y$n_patients,
    n_states = y$n_states,
    n_transitions = y$n_transitions,
    n_times = y$n_times
  ))
}

make_id_data_table <- function(x) {
  all_id_vars <- c("sample", "strategy_id", "patient_id", "state_id", 
                   "transition_id", "time_id")
  id_vars <- all_id_vars[all_id_vars %in% names(x)] # ID variables used
  
  id_dt <- as.data.table(x[id_vars])
  
  if ("time_id" %in% colnames(id_dt)) {
    ti <- x$time_intervals[match(id_dt$time_id, x$time_intervals$time_id)]
    ti <- ti[, c("time_start", "time_stop"), with = FALSE]
    id_dt <- cbind(id_dt, ti)
  }
  setattr(id_dt, "id_vars", id_vars[id_vars != "sample"])
  return(id_dt)
}

# Labels -----------------------------------------------------------------------
#' Set value labels
#' 
#' Update existing variables or create new ones that replace existing values
#' with more informative labels as in [`factor()`]. All modifications are performed 
#' by reference (see [`data.table::set()`] for more information about assignment by 
#' reference).
#' @param x A `data.table`.
#' @param labels A list of named vectors containing the values and labels of 
#' variables. The elements of each vector are the values of a variable and the 
#' names are the labels. The names of the list are the names of the variables.
#' See the output returned by [`get_labels()`] for an example.
#' @param new_names A character vector of the same length as `labels` where
#' each element denotes the name of a new variable to create for the
#' corresponding element in `labels`. If `NULL`, then the variables in `labels`
#' are modified and no new variables are created; otherwise, the existing variables
#' are not modified and new variables are created instead.
#' @param as_factor If `TRUE` factor variables are created; otherwise character
#' vectors are created. 
#' @return `x` is modified by reference and returned invisibly.  
#' @examples 
#' library("data.table")
#' labs <- list("strategy_id" = c("s1" = 1, 
#'                                "s2" = 2),
#'             "grp_id" = c("g1" = 1, 
#'                          "g2" = 2))
#' d1 <- data.table(strategy_id = 1:2, grp_id = 1:2)
#' d2 <- copy(d1); d3 <- copy(d2)
#' set_labels(d2, labels = labs)
#' set_labels(d3, labels = labs, new_names = c("strategy_name", "grp_name"))
#' d1
#' d2
#' d3
#' @seealso [`get_labels()`]
#' @export
set_labels <- function(x, labels, new_names = NULL, as_factor = TRUE) {
  
  labels <- labels[names(labels) %in% colnames(x)]
  
  if (length(labels) > 0) {
    for (i in 1:length(labels)) {
      old_name <- names(labels)[i]
      if (!is.null(new_names)) new_name <- new_names[i] else new_name <- old_name
      if (is.null(names(labels[[i]]))) {
        stop("Each element of 'labels' must be a named vector.")
      }
      x[,  (new_name) := factor(x[[old_name]], 
                                levels = labels[[i]],
                                labels = names(labels[[i]]))]
      if (!as_factor) x[, (new_name) := as.character(x[[new_name]])]
    }
  }
  invisible(x[])
}

#' Get value labels
#' 
#' Get value labels for the ID variables in a `hesim_data` object and create a list
#' of named vectors that can be passed to formatting and plotting functions. This
#' lets users create nice labels for treatment strategies, subgroups, health states,
#' and/or transitions when presenting results. 
#' @param object An object of class `hesim_data` created with [`hesim_data()`].
#' @param strategy The name of the column in the `strategy` element of `object`
#' containing labels for `strategy_id`.
#' @param grp The name of the column in the `patient` element of `object`
#' containing labels for `grp_id`.
#' @param state The name of the column in the `state` element of `object`
#' containing labels for `state_id`.
#' @param transition The name of the column in the `transition` element of `object`
#' containing labels for `transition_id`.
#' @param death_label The label to use for the death health state. By default a
#' label named "Death" will be concatenated to the labels for the non-death health
#' states. The death state can be omitted from labels for the health states by setting
#' `death_label = NULL`.
#' @return A list of named vectors containing the values and labels of 
#' variables. The elements of each vector are the values of a variable and the names 
#' are the labels. The names of the list are the names of the ID variables. 
#' @examples
#' library("data.table")
#' strategies <- data.table(
#'   strategy_id = c(1, 2),
#'   strategy_name = c("Strategy 1", "Strategy 2")
#' )
#' patients <- data.table(
#'   patient_id = seq(1, 4),
#'   age = c(50, 55, 60, 65),
#'   grp_id = c(1, 1, 2, 2),
#'   grp_name = rep(c("Age 50-59", "Age 60-69"), each = 2)
#' )
#' states <- data.table(
#'   state_id =  seq(1, 2),
#'   state_name = c("State 1", "State 2")
#' )
#' hesim_dat <- hesim_data(
#'   strategies = strategies,
#'   patients = patients,
#'   states = states
#' )
#' labs <- get_labels(hesim_dat)
#' labs
#' 
#' # Pass to set_labels()
#' d <- data.table(strategy_id = c(1, 1, 2, 2),
#'                 grp_id = c(1, 2, 1, 2))
#' set_labels(d, labs, new_name = c("strategy_name", "grp_name"))
#' d
#' @seealso [`hesim_data()`], [`set_labels()`]
#' @export
get_labels <- function(object, strategy = "strategy_name",
                       grp = "grp_name", state = "state_name",
                       transition = "transition_name", 
                       death_label = "Death") {
  check_is_class(object, "hesim_data", "object")

  # All possible ID variables
  tables <- c("strategies", "patients", "states", "transitions")
  id_vars <- c("strategy_id", "grp_id", "state_id", "transition_id")
  label_vars <- list(strategy, grp, state, transition)
  
  # Remove NULL labels and tables
  label_keep1 <- which(sapply(label_vars, function (z) !is.null(z)))
  table_keep <- which(tables %in% names(object))
  keep <- intersect(label_keep1, table_keep)
  if (length(keep) == 0) stop("There are no labels to get.")
  
  # Create table of non-NULL labels and tables
  m <- data.table(
    table = tables[keep], 
    id = id_vars[keep],
    label = unlist(label_vars[keep])
  )

  # Create labels
  create_labels <- function(object, id_var, label_var, table_name) {
    if (label_var %in% colnames(object[[table_name]])) {
      x <- as.data.table(object[[table_name]])
      x <- unique(x[, c(id_var, label_var), with = FALSE])
      if (length(unique(x[[id_var]])) != nrow(x)) {
        stop("There should be exactly one label for each ID value.",
             call. = FALSE)
      }
      v <- x[[id_var]]
      names(v) <- x[[label_var]]
      return(v)
    } else{
      return(NULL)
    } 
  }
  
  l <- vector(mode = "list", length = nrow(m))
  names(l) <- m$id
  for (i in 1:length(l)){
    labs <- create_labels(object, id_var = m$id[i], label_var = m$label[i],
                           table_name = m$table[i])
    if (!is.null(labs)) l[[i]] <- labs
  }
  l <- l[lengths(l) != 0]
  if (length(l) == 0) stop("The selected labels are not contained in the tables of 'object'.")
  
  # Add death label
  if ("state_id" %in% names(l) & !is.null(death_label)) {
    l$state_id <- c(l$state_id, max(l$state_id) + 1L)
    names(l$state_id)[length(l$state_id)] <- death_label
  }

  # Return
  return(l)
}