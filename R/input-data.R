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
#' \item{treatment_id}{Treatment id for treatment used at a given line of therapy within a treatment strategy.}
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
#' \item{transition_id}{Health state transition id.}
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
#' Each table must either be a \code{\link{data.frame}} or \code{\link{data.table}}. All id variables within 
#' each table must be numeric vectors of integers. 
#' @param strategies A table of treatment strategies. 
#' Must contain the column \code{strategy_id} denoting a unique strategy. Other columns are variables
#'  describing the characteristics of a treatment strategy. 
#' @param patients A table of patient observations. 
#' Must contain the column \code{patient_id} denoting a unique patient. The 
#' number of rows should be equal to the number of patients in the model.
#' Other columns are variables describing the characteristics of a patient.
#' @param lines A table of treatment lines used for each treatment strategy. Must contain the columns
#' \code{strategy_id}, denoting a treatment strategy, and \code{line}, denoting a treatment line. Other 
#' columns are variables describing the characteristics of a treatmnet line for a given treatment
#' strategy. A column denoting the treatment used for a given strategy and line would often
#' be specified.
#' @param states A table of health states. Must contain the column
#' \code{state_id}, which denotes a unique health state. The number of rows should
#' be equal to the number of health states in the model. Other columns can describe the
#' characteristics of a health state.
#' @param transitions A table of health state transitions. Must contain the column
#' \code{transition_id}, which denotes a unique transition; \code{from}, which denotes
#' the starting health state; and \code{to} which denotes the state that will be
#' transitioned to.
#' @return Returns an object of class "hesim_data", which is a list of data tables for
#' health economic simulation modeling.
#' @examples 
#' dt_strategies <- data.frame(strategy_id = c(1, 2))
#' dt_patients <- data.frame(patient_id = seq(1, 3), age = c(65, 50, 75),
#'                           gender = c("Female", "Female", "Male"))
#' dt_lines <- create_lines_dt(list(c(1, 2, 5), c(1, 2)))
#' dt_states <- data.frame(state_id =  seq(1, 3),
#'                         state_var = c(2, 1, 9))
#' hesim_dat <- hesim_data(strategies = dt_strategies,
#'                          patients = dt_patients,
#'                          states = dt_states,
#'                          lines = dt_lines)
#' @export
hesim_data <- function(strategies, patients, lines = NULL, states = NULL,
                          transitions = NULL){
  object <- new_hesim_data(strategies = strategies, patients = patients,
                           lines = lines, states = states, transitions = transitions)
  return(check(object))
}

new_hesim_data <- function(strategies, patients, lines = NULL, states = NULL,
                          transitions = NULL){
  data <- list()
  data$strategies <- strategies 
  data$patients <- patients
  data$lines <- lines
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
  
  # Lines
  if (!is.null(x$lines)){
      check_hesim_data_type(x$lines, "lines")
      if (!"strategy_id" %in% colnames(x$lines)){
        stop("'lines' must contain the column 'strategy_id'.", 
             call. = FALSE)
      }
      if (!"line" %in% colnames(x$lines)){
        stop("'lines' must contain the column 'line'.", 
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

#' Expand hesim_data
#' 
#' Expand the data tables from an object of class \code{\link{hesim_data}} into a single \code{\link{data.table}} 
#' of class "expanded_hesim_data". See the description of the return value for details on how this is done.
#' @param data An object of class \code{\link{hesim_data}}.
#' @param by A character vector of the names of the data tables in \code{\link{hesim_data}} to expand by.
#' @return This function is similar to \code{\link{expand.grid}}, but works for data frames or data tables. 
#' Specifically, it creates a \code{data.table} from all combinations of the supplied tables in \code{data}. 
#' The supplied tables are determined using the \code{by} argument. The resulting dataset is sorted by 
#' prioritizing id variables as follows: (i) \code{strategy_id}, (ii) \code{line}, (iii) \code{patient_id},
#' (iv) the health-related id variable (either \code{state_id} or \code{transition_id}).
#' @examples 
#' dt_strategies <- data.frame(strategy_id = c(1, 2))
#' dt_patients <- data.frame(patient_id = seq(1, 3), age = c(65, 50, 75),
#'                           gender = c("Female", "Female", "Male"))
#' dt_lines <- create_lines_dt(list(c(1, 2, 5), c(1, 2)))
#' dt_states <- data.frame(state_id =  seq(1, 3),
#'                         state_var = c(2, 1, 9))
#' hesim_dat <- hesim_data(strategies = dt_strategies,
#'                         patients = dt_patients,
#'                         states = dt_states,
#'                         lines = dt_lines)
#' expand_hesim_data(hesim_dat, by = c("strategies", "patients"))
#' @export
expand_hesim_data <- function(data, by = c("strategies", "patients")){
  if ("transitions" %in% by & "states" %in% by){
    stop("Cannot expand by both 'transitions' and 'states'.", call. = FALSE)
  }
  if (!all(by %in% names(hesim_data_sorting_map()))){
    stop("One of the elements specified in 'by' is not a table in 'hesim_data'.",
         call. = FALSE)
  }
  sorted_by <- hesim_data_sorted_by(by)
  tbl_list <- data[sorted_by]
  for (i in 1:length(tbl_list)){
    if (is.null(tbl_list[[i]])){
      stop("Cannot merge a NULL data table.")
    }
  }
  tbl_list <- lapply(tbl_list, function(x) as.data.frame(x))
  dat <- if ("lines" %in% sorted_by){
    tbl_list1 <- tbl_list[which(names(tbl_list) != "lines")]
    dat <- Reduce(function(...) merge(..., by = NULL), tbl_list1)
    dat <- data.table(merge(dat, tbl_list[[which(names(tbl_list) == "lines")]],
                     by = "strategy_id", sort = FALSE))
  } else{
    dat <- data.table(Reduce(function(...) merge(..., by = NULL, sort = FALSE),
                             tbl_list))
  }
  sort_hesim_data(dat, sorted_by)
  id_cols <- unlist(hesim_data_sorting_map()[sorted_by])
  nonid_cols <- colnames(dat)[!colnames(dat) %in% id_cols]
  dat <- dat[, c(id_cols, nonid_cols), with = FALSE]
  res <- list(data = dat,
              id_vars = unname(id_cols))
  class(res) <- "expanded_hesim_data"
  return(res)
}

hesim_data_sorting_map <- function(){
  list(strategies = "strategy_id",
       lines = "line",
       patients = "patient_id",
       transitions = "transition_id",
       states = "state_id")
}

hesim_data_sorted_by <- function(by){
  sorted_map <- hesim_data_sorting_map()
  sorted_by <- names(sorted_map)[names(sorted_map) %in% by]
  return(sorted_by)
}

sort_hesim_data <- function(data, sorted_by){
  setorderv(data, unlist(hesim_data_sorting_map()[sorted_by]))
}

# input_data class -------------------------------------------------------------
#' Input data for a statistical model
#' 
#' Create an object of class "input_data", which contains data for use
#' as an input to a statistical model. Data consists of (i) input matrices, \code{X},
#' (ii) id variables indexing the rows of each matrix in \code{X}, and (iii) the dimensions of
#' the \code{X} matrices. More details are provided under "Details" below. Note that an "input_data" 
#' object can often be created more easily using the generic function 
#' \code{\link{create_input_data}}. 
#' 
#' @param X A list of input matrices for predicting the values of each parameter in a statistical model. May also be
#' a list of lists of input matrices when a list of separate models is fit (e.g., with \link{flexsurvreg_list}).
#' @param strategy_id A numeric vector of integers denoting the treatment strategy represented by each row
#' in \code{X}.
#' @param n_strategies A scalar denoting the number of unique treatment strategies.
#' @param patient_id A numeric vector of integers denoting the patient represented by each row
#' in \code{X}.
#' @param n_patients A scalar denoting the number of unique patients.
#' @param line A numeric vector of integers denoting the treatment line represented by each row
#' in \code{X}. Not supported by currently available models.
#' @param n_lines A \code{\link{data.table}} denoting the number of treatment lines associated
#' with each treatment strategy. Should contain a column, "strategy_id", and a column,
#' "N". Not supported by currently available models.
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
#' @param time_fun A pointer to a C++ functor that can be used to update \code{X} as a function
#' of time in a simulation model. Not currently supported.
#' 
#' @details Each row of each matrix \code{X} is an input vector, \eqn{x_{hijk}}, where \eqn{h} denotes
#' a health-related index, \eqn{i} indexes a patient, \eqn{j} indexes a treatment line,
#' and \eqn{k} is a treatment strategy. A health-related index is either a health state
#' (e.g., \code{state_id}) or a transition between health states (e.g., \code{transition_id}).
#' In some cases, the health-related index \eqn{h} can be suppressed and separate models
#' can be fit for each health index. This is, for instance, the case in a partitioned survival 
#' model where separate models are fit for each survival endpoint. Likewise, models 
#' can be fit without multiple treatment lines as would, again, be the case in a 
#' partitioned survival analysis where sequential treatment would be incorporated by 
#' adding additional health states rather than by using the index \eqn{j}.
#' 
#' The rows of the matrices in \code{X} must be sorted in a manner consistent with the id variables.
#' The sorting order should be the same as specified in \code{\link{expand_hesim_data}}; that is,
#' the rows of \code{X} must be sorted by: (i) \code{strategy_id}, (ii) \code{line}, 
#' (iii) \code{patient_id}, and (iv) the health-related id variable (either \code{state_id} or
#'  \code{transition_id}).
#' @examples 
#' dt_strategies <- data.frame(strategy_id = c(1, 2))
#' dt_patients <- data.frame(patient_id = seq(1, 3), 
#'                           age = c(45, 47, 60),
#'                           female = c(1, 0, 0),
#'                           group = factor(c("Good", "Medium", "Poor")))
#' hesim_dat <- hesim_data(strategies = dt_strategies,
#'                         patients = dt_patients)
#' 
#' dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients"))$data
#' input_dat <- input_data(X = list(mu = model.matrix(~ age, dat)),
#'                        strategy_id = dat$strategy_id,
#'                        n_strategies = length(unique(dat$strategy_id)),
#'                        patient_id = dat$patient_id,
#'                        n_patients = length(unique(dat$patient_id)))
#' print(input_dat)
#' @export
input_data <- function(X, strategy_id, n_strategies,
                       patient_id, n_patients,
                       line = NULL, n_lines = NULL,
                       state_id = NULL, n_states = NULL,
                       transition_id = NULL, n_transitions = NULL,
                       time_fun = NULL){
  object <- new_input_data(X, strategy_id, n_strategies,
                                    patient_id, n_patients,
                                    line, n_lines,
                                    state_id, n_states,
                                    transition_id, n_transitions,
                                    time_fun)
  check(object)
  return(object)
}

new_input_data <- function(X, strategy_id, n_strategies,
                       patient_id, n_patients,
                       line = NULL, n_lines = NULL,
                       state_id = NULL, n_states = NULL,
                       transition_id = NULL, n_transitions = NULL,
                       time_fun = NULL){
  stopifnot(is.matrix(X) | is.list(X) | is.null(X))
  stopifnot(is.numeric(strategy_id))
  stopifnot(is.numeric(n_strategies))
  stopifnot(is.numeric(patient_id))
  stopifnot(is.numeric(n_patients))
  stopifnot(is.numeric(line) | is.null(line))
  stopifnot(is.data.table(n_lines) | is.null(n_lines))
  stopifnot(is.numeric(state_id) | is.null(state_id))
  stopifnot(is.numeric(n_states) | is.null(n_states))
  stopifnot(is.numeric(transition_id) | is.null(n_transitions))
  stopifnot(is.numeric(n_transitions) | is.null(n_transitions))
  if(!is.null(time_fun)){
   if(!inherits(time_fun, "externalptr")){
     stop("If not NULL, 'time_fun' must be of class 'externalptr'.",
          call. = FALSE)
   } 
  }
  
  object <- list(X = X,
                 strategy_id = strategy_id, n_strategies = n_strategies,
                 patient_id = patient_id, n_patients = n_patients,
                 line = line, n_lines = n_lines,
                 state_id = state_id, n_states = n_states,
                 transition_id = transition_id, n_transitions = n_transitions,
                time_fun = time_fun)
  object[sapply(object, is.null)] <- NULL
  class(object) <- "input_data"
  return(object)
}

#' @rdname check
check.input_data <- function(object){
  if (!is.list(object$X)){
    stop("'X' must be a list or a list of lists.", call. = FALSE)
  }
  X <- flatten_lists(object$X)
  if(!all(sapply(X, function(y) inherits(y, "matrix")))){
    stop("'X' must be a list or list of lists of matrices.", call. = FALSE)
  }
  X_nrows <- sapply(X, nrow)
  if (is.list(object$X)){
    if (!all(X_nrows[1] == X_nrows)){
      stop("The number of rows in each matrix in 'X' must be the same.",
           call. = FALSE)
    }
  }
  
  # Check that strategies in n_lines are correct
  if (!is.null(object$n_lines)){
    if(!all(sort(object$n_lines$strategy_id) == sort(unique(object$strategy_id)))){
      msg <- "'strategy_id' in 'n_lines' is not consistent with 'strategy_id'."
      stop(msg, call. = FALSE)
    }
  }
  
  id_vars <- c("strategy_id", "patient_id", "line", "state_id", "transition_id")
  id_vars_n <- c("n_strategies", "n_patients", "n_lines", "n_states", "n_transitions")
  for (i in 1:length(id_vars)){
    if (!is.null(object[[id_vars[i]]])){
      
      ## Check number of rows
      if(length(object[[id_vars[i]]]) != X_nrows[1]){
        msg <- paste0("The length of '", id_vars[i], "' does not equal the number of rows in the ",
                      "'X' matrices.")
        stop(msg, call. = FALSE)
      }
      
      # Check that n_strategies, n_patients, ..., is correct
        if (id_vars[i] != "line"){
          if(length(unique(object[[id_vars[i]]])) != object[id_vars_n[i]]){
            msg <- paste0("The number of unique observations in '", id_vars[i], 
                      "' does not equal '", id_vars_n[i], "'.")
            stop(msg, call. = FALSE)
          } 
          } else{
            line <- object[[id_vars[i]]]
            expected_n <- data.table("strategy_id" = object[["strategy_id"]],
                                     "line" = line)
            expected_n <- expected_n[, list(N = length(unique(line))), by = "strategy_id"]
            if(!all(expected_n$N == object[[id_vars_n[i]]]$N)){
              msg <- paste0("The number of unique observations in 'line'", 
                        "' within each 'strategy_id' does not equal 'n_lines$N'")
              stop(msg, call. = FALSE)
            }  
        } # end if else for line id variable
    } # end loop of id.vars
  }
  
  # Check if id variables are sorted properly 
  indices_df <- data.table(do.call("cbind", object[id_vars]))
  sorted_seq <- seq_len(nrow(indices_df))
  indices_df[, "row_num" := sorted_seq]
  by <- id_vars[sapply(object[id_vars], function(x) !is.null(x))]
  sort_hesim_data(indices_df, sorted_by = hesim_data_sorted_by(by))
  if(!all(indices_df$row_num == sorted_seq)){
    msg <- paste0("The id variables are not sorted correctly. The sort priority of the ",
                  "id variables must be as follows: 'strategy_id', 'line', 'patient_id' and ",
                   "the health-related id variable ('state_id' or 'transition_id').")
    stop(msg, call. = FALSE)
  }

  # Check if the number of unique observations is correct within groups
  indices_df[, "row_num" := NULL]
  for (i in 2:ncol(indices_df)){
    dt_by <- colnames(indices_df)[i - 1]
    col <- colnames(indices_df)[i]
    if (id_vars_n[i] == "n_lines"){
      len <- indices_df[, list(len = length(unique(get(col)))), 
                        by = c("strategy_id", dt_by)]$len
    } else{
      len <- indices_df[, list(len = length(unique(get(col)))), by = dt_by]$len 
    }
    if (id_vars_n[i] == "n_lines"){
      user_n <- object[[id_vars_n[i]]]$N
    } else{
      user_n <- object[[id_vars_n[i]]]
    }
    if (!all(unique(len) == user_n)){
      if (id_vars_n[i] == "n_lines"){
        msg <- paste0("The number of unique '", col, "' observations within each value",
                    " of 'strategy_id' and '", dt_by, " ' must be consistent with '", id_vars_n[i], "'.")
        stop(msg, call. = FALSE)
      } else{
        msg <- paste0("The number of unique '", col, "' observations within each value",
                    " of '", dt_by, " ' must equal '", id_vars_n[i], "'.")
        stop(msg, call. = FALSE)
      }
    }
  }
}

# Create input data from a fitted model ----------------------------------------
size_id_map <- function(){
  c(strategy_id = "n_strategies", 
    patient_id = "n_patients",
    line = "n_lines",
    state_id = "n_states",
    transition_id = "n_transitions")
}

get_input_data_id_vars <- function(data){
  map <- size_id_map()
  res <- list() 
  id_vars <- data$id_vars
  for (i in 1:length(id_vars)){
    res[[id_vars[i]]] <- data$data[[id_vars[i]]]
    if (id_vars[i] != "line"){
       res[[map[id_vars[i]]]] <- length(unique(data$data[[id_vars[i]]]))
    } else{
      n_lines <- data$data[, .N, by = c("strategy_id", "line")][, .N, by = "strategy_id"]
      res[[map[id_vars[i]]]] <- n_lines
    }
  }
  return(res)
}

#' Check data argument for \code{create_input_data} 
#' 
#' Check that data argument for \code{create_input_data} exists and that it is
#' of the correct type. 
#' @param data An object of class "expanded_hesim_data" returned by the function
#'  \code{\link{expand_hesim_data}}. 
#' @return If all tests passed, returns nothing; otherwise, throws an exception.
check_edata <- function(data){
  if (missing(data)){
    stop("'data' is missing with no default.")
  }
  if (!inherits(data, "expanded_hesim_data")){
    stop("'data' must be of class 'expanded_hesim_data'.")
  }   
}

#' Create input data
#' 
#' \code{create_input_data} is a generic function for creating an object of class
#' \code{\link{input_data}}. Model matrices are typically constructed based on the 
#' variables specified in the model \code{object} and the data specified in \code{data}, 
#' although there are some cases in which \code{\link{input_data}} can be created
#' from \code{object} alone.
#' @param object An object of the appropriate class. Currently supports
#' \code{\link{formula_list}}, \code{\link{lm}}, \code{\link{flexsurvreg}}, 
#'  \code{\link{flexsurvreg_list}}, and \code{\link{partsurvfit}}.
#' @param data An object of class "expanded_hesim_data" returned by the function
#'  \code{\link{expand_hesim_data}}. Used to look for the input variables needed to create an input matrix
#'  for use in a statistical models and the id variables for indexing rows in the input matrix. 
#' @param ... Further arguments passed to \code{\link{model.matrix}}.
#' @return An object of class \code{\link{input_data}}.
#' @examples 
#' library("flexsurv")
#' 
#' dt_strategies <- data.frame(strategy_id = c(1, 2))
#' dt_patients <- data.frame(patient_id = seq(1, 3), 
#'                           age = c(45, 47, 60),
#'                           female = c(1, 0, 0),
#'                           group = factor(c("Good", "Medium", "Poor")))
#' dt_states <- data.frame(state_id =  seq(1, 3),
#'                         state_name = factor(paste0("state", seq(1, 3))))
#' hesim_dat <- hesim_data(strategies = dt_strategies,
#'                         patients = dt_patients,
#'                         states = dt_states)
#'
#' # Class "lm"
#' expanded_dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "states"))
#' fit_lm <- stats::lm(costs ~ female + state_name, psm4_exdata$costs$medical)
#' input_dat <- create_input_data(fit_lm, expanded_dat)
#' class(input_dat)
#'
#' # Class "flexsurvreg"
#' expanded_dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients"))
#' fit_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "weibull")
#' input_dat <- create_input_data(fit_wei, expanded_dat)
#' class(input_dat)
#' @export
#' @rdname create_input_data
create_input_data <- function (object, ...) {
  if (missing(object)){
    stop("'object' is missing with no default.")
  }
  UseMethod("create_input_data", object)
}

formula_list_rec <- function(object, data, ...){
  x <- vector(mode = "list", length = length(object))
  names(x) <- names(object)
  for (i in 1:length(x)){
    if (inherits(object[[i]], "formula")){
      x[[i]] <- stats::model.matrix(object[[i]], data = data$data, ...)
    } else{
      x[[i]] <- formula_list_rec(object[[i]], data = data, ...)
    }
  }
  return(x)
}

#' @export
#' @rdname create_input_data
create_input_data.formula_list <- function(object, data, ...){
  check_edata(data)
  X_list <- formula_list_rec(object, data, ...)
  args <- c(list(X = X_list),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

get_terms <- function(object){
  tt <- stats::terms(object)
  return(stats::delete.response(tt))
}

#' @export
#' @rdname create_input_data
create_input_data.stateval_means <- function(object, ...){
  # Create expanded hesim data
  n_states <- dim(object$values)[2]
  states_dt <- data.table(state_id = seq(1, n_states))
  patients_dt <- data.table(patient_id = object$patient_id)
  strategies_dt <- data.table(strategy_id = object$strategy_id)
  hesim_dat <- hesim_data(strategies = strategies_dt,
                          patients = patients_dt,
                          states = states_dt)
  edat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "states"))   
  
  # Create input data
  args <- c(list(X = NULL),
           get_input_data_id_vars(edat))
  return(do.call("new_input_data", args))
}

#' @export 
#' @rdname create_input_data
create_input_data.lm <- function(object, data, ...){
  check_edata(data)
  terms <- get_terms(object)
  X <- stats::model.matrix(terms, data = data$data, ...)
  args <- c(list(X = list(mu = X)),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

#' @export 
#' @rdname create_input_data
create_input_data.lm_list <- function(object, data, ...){
  check_edata(data)
  X_list <- vector(mode = "list", length = length(object))
  names(X_list) <- names(object)
  for (i in 1:length(X_list)){
    terms <- get_terms(object[[i]])
    X_list[[i]] <- list(mu = stats::model.matrix(terms, data = data$data, ...))
  }
  args <- c(list(X = X_list),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

create_input_data_flexsurvreg_X <- function(object, data, ...){
  pars <- object$dlist$pars
  X_list <- vector(mode = "list", length = length(pars))
  names(X_list) <- pars
  for (i in 1:length(pars)){
    form <- object$all.formulae[[pars[i]]]
    if (is.null(form)){
      form <- stats::formula(~1)
    } else{
      form <- stats::delete.response(stats::terms(form))
    }
    X_list[[i]] <- stats::model.matrix(form, data = data$data, ...)
  }
  return(X_list)
}

#' @export
#' @rdname create_input_data
create_input_data.flexsurvreg <- function(object, data,...){
  check_edata(data)
  X_list <- create_input_data_flexsurvreg_X(object, data, ...)
  args <- c(list(X = X_list),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

#' @export
#' @rdname create_input_data
create_input_data.flexsurvreg_list <- function(object, data,...){
  check_edata(data)
  X_list_2d <- vector(mode = "list", length = length(object))
  names(X_list_2d) <- names(object)
  for (i in 1:length(object)){
    X_list_2d[[i]] <- create_input_data_flexsurvreg_X(object[[i]], data, ...)
  }
  args <- c(list(X = X_list_2d),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

#' @export
#' @rdname create_input_data
create_input_data.partsurvfit <- function(object, data, ...){
  check_edata(data)
  return(create_input_data.flexsurvreg_list(object$models, data, ...))
}

#' @export
create_input_data.joined_flexsurvreg_list <- function(object, data,...){
  check_edata(data)
  models <- object$models
  X_list_3d <- vector(mode = "list", length = length(models))
  names(X_list_3d) <- names(models)
  for (i in 1:length(models)){
    X_list_3d[[i]] <- vector(mode = "list", length = length(models[[i]]))
    names(X_list_3d[[i]]) <- names(models[[i]])
    for (j in 1:length(models[[i]])){
      X_list_3d[[i]][[j]] <- create_input_data_flexsurvreg_X(models[[i]][[j]], data, ...)
    }
  }
  args <- c(list(X = X_list_3d),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}
