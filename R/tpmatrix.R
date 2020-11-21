# Transition probability matrix ------------------------------------------------
transform_dots <- function(dots, data){
  n_dots <- length(dots)
  for (i in 1:n_dots){
    name_i <- names(dots)[i]
    data[[name_i]] <-  eval(dots[[i]], data)
  }
  return(data[names(dots)])
}

#' Names for elements of a transition probability matrix
#' 
#' Create names for all elements of a transition probability matrix given
#' names for the health states. This is useful for flattening a transition
#' probability matrix (rowwise) into a vector and naming the resulting vector. 
#' The name of an element of the flattened vector representing a transition from
#' the ith state to the jth state is of the form 
#' `paste0(prefix, state_i, sep, state_j)`.
#' 
#' @param states A character vector of the names of health states in the 
#' transition matrix.
#' @param prefix A prefix that precedes the described transitions between states.
#' @param sep A character string to separate the terms representing
#' state `i` and state `j`.
#' @examples 
#' tpmatrix_names(LETTERS[1:4])
#' tpmatrix_names(LETTERS[1:4], prefix = "")
#' tpmatrix_names(LETTERS[1:4], prefix = "", sep = ".")
#' 
#' @return A character vector containing a name for each element of the transition
#' probability matrix encompassing all possible transitions. 
#' 
#' @export
tpmatrix_names <- function(states, prefix = "p_", sep = "_"){
  x <- paste0(prefix, states)
  y <- paste0(sep, states)
  return(c(t(outer(x, y, paste0))))
}

define_tpmatrix <- function(...){
  x <- as.list(substitute(c(...))[-1])
  n_states <- sqrt(length(x))
  if (!is_whole_number(n_states)){
    stop("tpmatrix() must be a square matrix.", call. = FALSE)
  } 
  names(x) <- tpmatrix_names(states = paste0("s", 1:n_states),
                             prefix = "")
  which_C <- which(x == "C")
  x[which_C] <- NA
  attr(x, "complement") <- rep(0, length(x))
  attr(x, "complement")[which_C] <- 1
  class(x) <- "tpmatrix_expr"
  return(x)
}

get_matind <- function(n_cols, i, j){
  ind <- (i - 1) * n_cols + j
  return(as.integer(ind))
}

get_matrow <- function(x, i, n_states){
  start <- get_matind(n_states, i, 1)
  return(x[, start:(start + n_states - 1)])
}

replace_C <- function(x, complement){
  
  n_states <- sqrt(ncol(x))
  for (i in 1:n_states){
    x_i <- get_matrow(x, i, n_states)
    complement_i <- complement[get_matind(n_states, i, 1:n_states)]
    if (sum(complement_i) > 1){
      stop("Only one 'C' is allowed per matrix row.", call. = FALSE)
    }
    if (sum(complement_i) == 1){
      C_val <- 1 - rowSums(x_i, na.rm = TRUE)
      C_index <- which(complement_i == 1)
      set(x, j = get_matind(n_states, i, C_index), value = C_val)
    }
  }
}

replace_Qdiag <- function(x, n_states) {
  for (i in 1:n_states){
    x_i <- get_matrow(x, i, n_states)
    q_diag_i <- -rowSums(x_i, na.rm = TRUE)
    q_diag_index_i <- get_matind(n_states, i, i)
    x[, q_diag_index_i] <- q_diag_i
  }
  return(x)
}

#' Transition probability matrix
#' 
#' `tpmatrix()` both defines and evaluates a transition probability matrix in which 
#' elements are expressions. It can be used within [define_tparams()] to 
#' create a transition probability matrix or directly to create a [tparams_transprobs()] 
#' object. These are, in turn, ultimately used to create a [CohortDtstmTrans] object
#' for simulating health state transitions.
#' 
#' @param ... Named values of expressions defining elements of the matrix. See
#' "Details" and the example below.
#' 
#' @details The matrix is filled rowwise, meaning that each row should sum to 1.
#' The complementary probability (equal to 1 minus the sum of the probabilities
#'  of all other rows) can be conveniently referred to as `C`. 
#' 
#' @return Returns a `tpmatrix` object that inherits from `data.table`
#' where each column is an element of the transition probability matrix with
#' elements ordered rowwise. 
#' 
#' @examples 
#' p <- c(.7, .6)
#' tpmatrix(
#'   C, p,
#'   0, 1
#' )
#' @seealso [define_model()], [define_tparams()], 
#' [tpmatrix_id()], [tparams_transprobs()], [CohortDtstmTrans()]
#' @export
tpmatrix <- function(...){
  m_def <- define_tpmatrix(...)
  m <- as.data.table(transform_dots(m_def, as.list(parent.frame())))
  replace_C(m, attr(m_def, "complement"))
  setattr(m, "class", c("tpmatrix", "data.table", "data.frame"))
  return(m)
}

#' Transition probability matrix IDs
#' 
#' Creates ID variables for each row returned by `tpmatrix()`. This function is
#'  most conveniently used along with `tpmatrix()` to construct a 
#'  `tparams_transprobs()` object.
#' 
#' 
#' @param object An object of class `expanded_hesim_data` returned by 
#' `expand.hesim_data()`. This dataset must be expanded by treatment
#' strategies, patients, and optionally time intervals.
#' @param n_samples The number of parameters samples used for the probabilistic
#' sensitivity analysis (PSA).
#' 
#' @return Returns a `data.table` with the same columns in `object` repeated
#' `n_samples` times. That is, to facilitate creation of a `tparams_transprobs()`
#' object,  there is one row for each parameter sample,
#' treatment strategy, patient, and optionally time interval.
#' 
#' @examples 
#' strategies <- data.frame(strategy_id = c(1, 2))
#' patients <- data.frame(patient_id = seq(1, 3), age = c(65, 50, 75),
#'                         gender = c("Female", "Female", "Male"))
#' hesim_dat <- hesim_data(strategies = strategies,
#'                         patients = patients)
#' input_data <- expand(hesim_dat, by = c("strategies", "patients"))    
#' tpmatrix_id(input_data, n_samples = 2)                   
#' @seealso [tpmatrix()], [tparams_transprobs()], [expand.hesim_data()]
#' @export
tpmatrix_id <- function(object, n_samples){
  # Checks
  check_is_class(object, "expanded_hesim_data", "object")
  if (!identical(attr(object, "id_vars"), c("strategy_id", "patient_id")) &&
      !identical(attr(object, "id_vars"), c("strategy_id", "patient_id", "time_id"))
  ) {
    stop(paste0("'object' must be either be expanded by 'strategy_id', 'patient_id', ", 
                "and optionally 'time_id'."), 
         call. = FALSE)
  }
  
  # Create data.table
  id_names <- c("strategy_id", "patient_id", "grp_id", "patient_wt", 
                "time_id", "time_start", "time_stop")
  cols_to_keep <- match(id_names, colnames(object))
  cols_to_keep <- cols_to_keep[!is.na(cols_to_keep)]
  
  id_dt <- object[rep(1:nrow(object), times = n_samples), cols_to_keep, 
                  with = FALSE]
  id_dt[, sample := rep(1:n_samples, each = nrow(object))]
  setcolorder(id_dt, c("sample", colnames(id_dt)[!ncol(id_dt)]))
  setattr(id_dt, "class", c("tpmatrix_id", "data.table", "data.frame"))
  return(id_dt[, ])
}

# Transition intensity matrix --------------------------------------------------
#' Transition intensity matrix
#' 
#' `qpmatrix()` creates transition intensity matrices where elements represent
#' the instantaneous risk of moving between health states. 
#' 
#' @param q A two-dimensional object that can be passed to [as.matrix()] containing
#' elements of the transition intensity matrix. A column represents a transition
#' from state \eqn{r} to state \eqn{s}. Each row represents elements of a different
#' transition intensity matrix. See "Details" for more information.
#' 
#' @param trans_mat Just as in [IndivCtstmTrans], a transition matrix 
#' describing the states and transitions in a multi-state model.
#' 
#' @details The object `q` must only contain non-zero and non-diagonal elements
#' of a transition intensity matrix. The diagonal elements are automatically computed
#' as the negative sum of the other rows. 
#' 
#' @return Returns a `qmatrix` object where each row is a flattened transition 
#' intensity matrix with elements ordered rowwise.  
#' 
#' @examples 
#' # 3 state irreversible model
#' tmat <- rbind(c(NA, 1, 2),
#'               c(NA, NA, 3),
#'               c(NA, NA, NA)) 
#' q12 <- c(.8, .7)
#' q13 <- c(.2, .3)
#' q23 <- c(1.1, 1.2)
#' q <- data.frame(q12, q13, q23)
#' qmat <- qmatrix(q, trans_mat = tmat)
#' print(qmat)
#' 
#' # Compute matrix exponential of each row
#' matrix_exp(qmat)
#' 
#' @seealso [tpmatrix()]
#' @export
qmatrix <- function(q, trans_mat){
  q <- as.matrix(q)
  trans <- c(t(trans_mat))
  n_states <- nrow(trans_mat)
  qmat <- matrix(0, nrow = nrow(q), ncol = length(trans))
  qmat[, which(!is.na(trans))] <- q
  colnames(qmat) <- tpmatrix_names(states = paste0("s", 1:n_states),
                                   prefix = "")
  qmat <- replace_Qdiag(qmat, n_states)
  class(qmat) <- "qmatrix"
  attr(qmat, "n_states") <- n_states
  return(qmat)
}

# Matrix exponential -----------------------------------------------------------
#' Matrix exponential
#' 
#' This is a wrapper around [msm::MatrixExp()] that computes the matrix
#' exponential of multiple square matrices. 
#' 
#' @param x A two-dimensional array like object where each row is a square
#' matrix ordered rowwise.
#' @param t An optional scaling factor for `x`. 
#' @param ... Arguments to pass to [msm::MatrixExp].
#' 
#' @details This function is most useful when exponentiating transition intensity
#' matrices to produce transition probability matrices. To create transition
#' probability matrices for discrete time state transition models with annual
#' cycles, set `t=1`. An array of matrices is returned which can be used
#' to create the `value` element of a [tparams_transprobs] object. See 
#' [qmatrix()] for an example. 
#' 
#' @return Returns an array of exponentiated matrices. 
#' 
#' @seealso [qmatrix()]
#' @export
matrix_exp <- function(x, t = 1, ...) {
  x <- as.matrix(x)
  n <- nrow(x)
  n_states <- sqrt(ncol(x))
  if (!is_whole_number(n_states)) stop("Each row of 'x' must be square matrix.")
  res <- array(NA, dim = c(n_states, n_states, n * length(t)))
  for (i in 1:n) {
    xmat_i <- matrix(x[i, ], byrow = TRUE, nrow = n_states)
    res[,, i] <- msm::MatrixExp(xmat_i, t = t, ...)
  }
  return(res)
}