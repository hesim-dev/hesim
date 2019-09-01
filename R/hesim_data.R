# hesim data -------------------------------------------------------------------

#' Create a data table of treatment lines
#' 
#' Convert a list of treatment lines for multiple treatment strategies to a \code{\link{data.table}}
#' suitable for use with \link{hesim_data}.
#' @param strategy_list A list where each element is a treatment strategy 
#' consisting of a vector of treatments. 
#' @param strategy_ids A numeric vector denoting the numeric id of each strategy
#' in \code{strategy_list}.
#' @return Returns a \code{\link{data.table}} in tidy format with three columns
#' \describe{
#' \item{strategy_id}{Treatment strategy ids.}
#' \item{line}{Line of therapy.}
#' \item{treatment_id}{Treatment ID for treatment used at a given line of therapy within a treatment strategy.}
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
#' the states and transitions in a multi-state model suitable for use with \link{hesim_data}. 
#' @param trans_mat A transition matrix in the format from the \link[mstate]{mstate} package. 
#' See \link{IndivCtstmTrans}.
#' @return Returns a \code{\link{data.table}} in tidy format with three columns
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

#' Data for health-economic simulation modeling
#' 
#' A list of tables required for health-economic simulation modeling.
#' Each table must either be a \code{\link{data.frame}} or \code{\link{data.table}}. All ID variables within 
#' each table must be numeric vectors of integers. 
#' @param strategies A table of treatment strategies. 
#' Must contain the column \code{strategy_id} denoting a unique strategy. Other columns are variables
#'  describing the characteristics of a treatment strategy. 
#' @param patients A table of patient observations. 
#' Must contain the column \code{patient_id} denoting a unique patient. The 
#' number of rows should be equal to the number of patients in the model.
#' Other columns are variables describing the characteristics of a patient.
#' @param states A table of health states. Must contain the column
#' \code{state_id}, which denotes a unique health state. The number of rows should
#' be equal to the number of health states in the model. Other columns can describe the
#' characteristics of a health state.
#' @param transitions A table of health state transitions. Must contain the column
#' \code{transition_id}, which denotes a unique transition; \code{from}, which denotes
#' the starting health state; and \code{to} which denotes the state that will be
#' transitioned to.
#' @param times A table of time intervals. Must contain the column \code{time_start}, 
#' which denotes the starting time of times intervals defined between \code{time_start}
#'  and \code{time_stop}. \code{hesim_data()} automatically creates the columns
#'  \code{time_stop} and \code{time_id} (a vector of integers sorted by \code{time_start}).
#'  Time intervals are closed on the left and open on the right and the final time
#'   interval is defined from \code{time_start} to infinity.
#' @return Returns an object of class "hesim_data", which is a list of data tables for
#' health economic simulation modeling.
#' @seealso \code{\link{expand.hesim_data}}
#' @examples 
#' strategies <- data.frame(strategy_id = c(1, 2))
#' patients <- data.frame(patient_id = seq(1, 3), age = c(65, 50, 75),
#'                           gender = c("Female", "Female", "Male"))
#' states <- data.frame(state_id =  seq(1, 3),
#'                         state_var = c(2, 1, 9))
#' times <- data.frame(time_id = c(1, 2, 3),
#'                     time_start = c(0, 4, 9),
#'                     time_stop = c(4, 9, Inf))
#' hesim_dat <- hesim_data(strategies = strategies,
#'                          patients = patients,
#'                          states = states,
#'                          times = times)
#' @export
hesim_data <- function(strategies, patients, states = NULL,
                       transitions = NULL, times = NULL){
  object <- new_hesim_data(strategies = strategies, patients = patients,
                           states = states, transitions = transitions,
                           times = times)
  return(check(object))
}

new_hesim_data <- function(strategies, patients, states = NULL,
                           transitions = NULL, times = NULL){
  data <- list()
  data$strategies <- strategies 
  data$patients <- patients
  data$states <- states
  data$transitions <- transitions
  data$times <- times
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
  
  # Times
  if (!is.null(x$times)){
    check_hesim_data_type(x$times, "times")
    
    ## time_start    
    if (!"time_start" %in% colnames(x$times)){
      stop("'times' must contain the column 'time_start'.",
           call. = FALSE)
    }
    
    ## time_stop
    time_stop <- shift(x$times$time_start, type = "lead")
    time_stop[is.na(time_stop)] <- Inf
    x$times$time_stop <- time_stop
    
    ## time_id
    x$times$time_id <- 1:nrow(x$times)
    x$times <- setcolorder(x$times, c("time_id", "time_start", "time_stop"))
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
#' \code{\link{hesim_data}} objects.
#' @export
#' @keywords internal
expand <- function(object, by){
  UseMethod("expand", object)
}

#' Expand hesim_data
#' 
#' Expand the data tables from an object of class \code{\link{hesim_data}} into a data table in 
#' long format (one row for each combination of observations as specified with the 
#' ID variables from the tables specified with the \code{by} argument).
#'  See "Details" for an explanation of how the expansion is done.
#' @param object An object of class \code{\link{hesim_data}}.
#' @param by A character vector of the names of the data tables in \code{\link{hesim_data}} to expand by.
#' @details This function is similar to \code{\link{expand.grid}}, but works for data frames or data tables. 
#' Specifically, it creates a \code{data.table} from all combinations of the supplied tables in \code{object}. 
#' The supplied tables are determined using the \code{by} argument. The resulting dataset is sorted by 
#' prioritizing ID variables as follows: (i) \code{strategy_id}, (ii) \code{patient_id},
#' (iii) the health-related ID variable (either \code{state_id} or \code{transition_id}), and
#' (iv) the time interval (i.e., \code{time_id}).
#' @return An object of class "expanded_hesim_data", which is a \code{data.table} with an "id_vars" 
#' attribute containing the names of the ID variables in the data table.
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
#' @export
expand.hesim_data <- function(object, by = c("strategies", "patients")){
  if ("transitions" %in% by & "states" %in% by){
    stop("Cannot expand by both 'transitions' and 'states'.", call. = FALSE)
  }
  if (!all(by %in% names(hesim_data_sorting_map()))){
    stop("One of the elements specified in 'by' is not a table in 'hesim_data'.",
         call. = FALSE)
  }
  sorted_by <- hesim_data_sorted_by(by)
  tbl_list <- object[sorted_by]
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
#' Attributes for ID variables
#' 
#' Stores metadata related to the ID variables used to index input data and prediction objects.
#' 
#' @param strategy_id A numeric vector of integers denoting the treatment strategy represented by each row
#' in \code{X}.
#' @param n_strategies A scalar denoting the number of unique treatment strategies.
#' @param patient_id A numeric vector of integers denoting the patient represented by each row
#' in \code{X}.
#' @param n_patients A scalar denoting the number of unique patients.
#' @param state_id A numeric vector of integers denoting the health state represented by each row
#' in \code{X}.
#' @param n_states A scalar denoting the number of unique health states.
#' @param transition_id A numeric vector denoting the 
#' health state transition represented by each row in \code{X}. This must only be specified when
#' estimating the health state transitions with a joint likelihood function. If independent
#' models are fit for each transition, then separate \code{X} matrices must be specified
#' for each transition. Note that this is not currently supported but
#' will be supported once \code{hesim} provides support for state transition modeling.
#' @param n_transitions A scalar denoting the number of unique transitions. 
#' Not supported by currently available models. 
#' @param time_id A numeric vector of integers denoting a unique time interval.
#' @param time_intervals A \code{data.table} denotes unique time intervals. Must 
#' contain the columns \code{time_id}, \code{time_start}, and \code{time_stop}.
#' \code{time_start} is the starting time of an interval and \code{time_stop} is
#' the stopping time of an interval. Time intervals are closed on the left and
#' open on the right, and in the final interval, \code{time_stop} is equal to 
#' infinity. 
#' @param n_times A scalar denoting the number of time intervals. Equal to the
#' number of rows in \code{time_intervals}.
#' @param time_reset If \code{TRUE}, then time intervals reset each time a patient enters a new health 
#'   state. This is relevant if, for example, costs vary over time within health states. 
#'   If \code{FALSE} then time intervals are based on time since the start of the simulation.
#'  
#' @details The sorting order should be the same as specified in \code{\link{expand.hesim_data}}; that is,
#' observations must be sorted by: (i) \code{strategy_id}, (ii) \code{line}, 
#' (iii) \code{patient_id}, and (iv) the health-related ID variable (either \code{state_id} or
#'  \code{transition_id}).
#'  @export
id_attributes <- function(strategy_id, n_strategies,
                          patient_id, n_patients,
                          state_id = NULL, n_states = NULL,
                          transition_id = NULL, n_transitions = NULL,
                          time_id = NULL, time_intervals = NULL, n_times = NULL,
                          time_reset = NULL,
                          time_fun = NULL){
  object <- new_id_attributes(strategy_id, n_strategies,
                              patient_id, n_patients,
                              state_id, n_states,
                              transition_id, n_transitions,
                              time_id, time_intervals, n_times,
                              time_reset,
                              time_fun)
  check(object)
  return(object)
}

new_id_attributes <- function(strategy_id, n_strategies,
                              patient_id, n_patients,
                              state_id = NULL, n_states = NULL,
                              transition_id = NULL, n_transitions = NULL,
                              time_id = NULL, time_intervals = NULL, n_times = NULL,
                              time_reset = NULL,
                              time_fun = NULL){
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
  stopifnot(is.logical(time_reset) | is.null(time_reset))
  
  object <- list(strategy_id = strategy_id, n_strategies = n_strategies,
                 patient_id = patient_id, n_patients = n_patients,
                 state_id = state_id, n_states = n_states,
                 transition_id = transition_id, n_transitions = n_transitions,
                 time_id = time_id, time_intervals = time_intervals, n_times = n_times,
                 time_reset = time_reset,
                 time_fun = time_fun)
  object[sapply(object, is.null)] <- NULL
  class(object) <- "id_attributes"
  return(object)
}

#' @rdname check
check.id_attributes <- function(object){
  id_vars <- c("strategy_id", "patient_id", "state_id", "transition_id",
               "time_id")
  id_vars_n <- c("n_strategies", "n_patients", "n_states", "n_transitions",
                 "n_times")
  for (i in 1:length(id_vars)){
    if (!is.null(object[[id_vars[i]]])){
      # Check that n_strategies, n_patients, ..., is correct
      if(length(unique(object[[id_vars[i]]])) != object[id_vars_n[i]]){
        msg <- paste0("The number of unique observations in '", id_vars[i], 
                      "' does not equal '", id_vars_n[i], "'.")
        stop(msg, call. = FALSE)
      } 
    } # end loop of id.vars
  }
  
  # Check if id variables are sorted properly 
  indices_df <- data.table(do.call("cbind", object[id_vars]))
  sorted_seq <- seq_len(nrow(indices_df))
  indices_df[, "row_num" := sorted_seq]
  by <- id_vars[sapply(object[id_vars], function(x) !is.null(x))]
  sort_hesim_data(indices_df, sorted_by = hesim_data_sorted_by(by))
  if(!all(indices_df$row_num == sorted_seq)){
    msg <- paste0("The ID variables are not sorted correctly. The sort priority of the ",
                  "ID variables must be as follows: 'strategy_id', 'patient_id' ",
                  "the health-related ID variable ('state_id' or 'transition_id') ",
                  "and 'time_id'.")
    stop(msg, call. = FALSE)
  }
  
  # Check if the number of unique observations is correct within groups
  indices_df[, "row_num" := NULL]
  for (i in 2:ncol(indices_df)){
    dt_by <- colnames(indices_df)[i - 1]
    col <- colnames(indices_df)[i]
    len <- indices_df[, list(len = length(unique(get(col)))), 
                      by = c("strategy_id", dt_by)]$len
    user_n <- object[[id_vars_n[i]]]
    if (!all(unique(len) == user_n)){
      msg <- paste0("The number of unique '", col, "' observations within each value",
                    " of '", dt_by, " ' must equal '", id_vars_n[i], "'.")
      stop(msg, call. = FALSE)
    }
  }
}