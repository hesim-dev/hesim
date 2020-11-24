# 3D array ---------------------------------------------------------------------
#' Convert between 2D tabular objects and 3D arrays
#' 
#' Convert a 2-dimensional tabular object where each row stores a flattened
#' square matrix to a 3-dimensional array of square matrices and vice versa. 
#' This allows multiple transition matrices to be stored as either tabular objects 
#' (e.g., matrices, data frames, etc) or as arrays. 
#' 
#' @param x For `as_array3()` a 2-dimensional tabular object where each row stores a flattened
#' square matrix ordered rowwise. Reasonable classes are `matrix`, `data.frame`,
#' `data.table`, and `tpmatrix`. For `as_tpmatrix()` a 3-dimensional array 
#' where each slice is a square matrix.
#' 
#'
#' @examples 
#' p_12 <- c(.7, .6)
#' pmat <- tpmatrix(
#'  C, p_12,
#'  0, 1
#' )
#' pmat
#' 
#' as_array3(pmat)
#' as_array3(as.matrix(pmat))
#' as_tpmatrix(as_array3(pmat))
#' @return For `as_array3()` a 3-dimensional array of square matrices; 
#' for `as_tpmatrix()` a [tpmatrix] object.
#' 
#' @seealso [tpmatrix]
#' @export
as_array3 <- function(x) {
  if (length(dim(x)) != 2) stop("'x' must be a 2-dimensional object.")
  n_states <- sqrt(ncol(x))
  y <- aperm(array(c(t(x)),
                   dim = c(n_states, n_states, nrow(x))),
             c(2, 1, 3))
  return(y)
}

#' @name as_array3
#' @export
as_tpmatrix <- function(x) {
  n_col <- dim(x)[1] * dim(x)[2]
  y <- matrix(c(aperm(x, c(2, 1, 3))), ncol = n_col, byrow = TRUE)
  colnames(y) <- tpmatrix_names(states = paste0("s", 1:dim(x)[1]),
                                prefix = "")
  return(tpmatrix(y))
}

# Transition probability matrix ------------------------------------------------
eval_dots <- function(dots, data){
  n_dots <- length(dots)
  res <- vector(mode = "list", length = n_dots)
  complement <- vector(mode = "list", length = n_dots)
  
  update_complement <- function(complement, x) {
    if (is.null(ncol(x))) {
      return(complement)
    } else{
      return(rep(complement, ncol(x)))
    }
  }
  
  for (i in 1:n_dots){
    res[[i]] <-  eval(dots[[i]], data)
    complement[[i]] <- update_complement(attr(dots, "complement")[i], res[[i]])
  }
  attr(res, "complement") <- unlist(complement)
  return(res)
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
#' @param ... Named values of expressions defining elements of the matrix. Each
#' element of `...` should either be a vector or a 2-dimensional tabular object 
#' such as a data frame. See "Details" and the examples below.
#' 
#' @details A `tpmatrix` is a 2-dimensional tabular object that stores flattened
#' square transition probability matrices in each row. Each transition probability
#' matrix is filled rowwise. The complementary probability 
#' (equal to \eqn{1} minus the sum of the probabilities
#'  of all other elements in a row of a transition probability matrix)
#'  can be conveniently referred to as `C`. There can 
#'  only be one complement for each row in a transition 
#'  probability matrix.
#' 
#' @return Returns a `tpmatrix` object that inherits from `data.table`
#' where each column is an element of the transition probability matrix with
#' elements ordered rowwise. 
#' 
#' @examples 
#' # Pass vectors
#' p_12 <- c(.7, .6)
#' tpmatrix(
#'   C, p_12,
#'   0, 1
#' )
#' 
#' tpmatrix(
#'   C, p_12,
#'   C, 1
#' )
#' 
#' # Pass matrix
#' pmat <- matrix(c(.5, .5, .3, .7), byrow = TRUE, ncol = 4)
#' tpmatrix(pmat)
#' 
#' # Pass vectors and data frames
#' p1 <- data.frame(
#'   p_12 = c(.7, .6), 
#'   p_13 = c(.1, .2)
#' )
#'
#' p2 <- data.frame(
#'   p_21 = 0,
#'   p_22 = c(.4, .45),
#'   p_23 = c(.6, .55)
#' )
#'
#' p3 <- data.frame(
#'   p_31 = c(0, 0),
#'   p_32 = c(0, 0),
#'   p_33 = c(1, 1)
#' )
#'
#' tpmatrix(
#'   C, p1,
#'   p2,
#'   p3
#' )
#'
#' 
#' @seealso [define_model()], [define_tparams()], 
#' [tpmatrix_id()], [tparams_transprobs()], [CohortDtstmTrans()]
#' @export
tpmatrix <- function(...){
  # Evaluate
  m_def <- define_tpmatrix(...)
  m <- eval_dots(m_def, as.list(parent.frame()))
  complement <- attr(m, "complement")
  m <- as.data.table(m)

  # Some checks
  n_states <- sqrt(ncol(m))
  if (!is_whole_number(n_states)){
    stop("tpmatrix() must be a square matrix.", call. = FALSE)
  }
  colnames(m) <- tpmatrix_names(states = paste0("s", 1:n_states),
                                prefix = "")

  # Replace complement
  replace_C(m, complement)
  
  # Return
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
#' @param q A two-dimensional tabular object that can be passed to [as.matrix()] containing
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
#' @return An array of transition intensity matrices with the third dimension 
#' equal to the number of rows in `q`.
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
#' # Matrix exponential of each matrix in array
#' expmat(qmat)
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
  return(as_array3(qmat))
}

# Matrix exponential -----------------------------------------------------------
#' Matrix exponential
#' 
#' This is a wrapper around [msm::MatrixExp()] that computes the exponential
#' of multiple square matrices. 
#' 
#' @param x An array of matrices. 
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
expmat <- function(x, t = 1, ...) {
  n <- dim(x)[3]
  n_states <- dim(x)[1]
  res <- array(NA, dim = c(n_states, n_states, n * length(t)))
  for (i in 1:n) {
    res[,, i] <- msm::MatrixExp(x[,, i], t = t, ...)
  }
  return(res)
}

# Relative risks ---------------------------------------------------------------
#' Apply relative risks to transition probability matrices
#' 
#' Elements of transition probability matrices are multiplied by relative risks
#' and the transition probability matrices are adjusted so that rows sum to 1. 
#' Operations are vectorized and each relative risk is multiplied by every 
#' transition matrix (stored in 3-dimensional arrays).  
#' 
#' @param x A 3-dimensional array where each slice is a square transition
#' probability matrix.
#' @param rr A 2-dimensional tabular object such as a matrix or data frame where each
#' column is a vector of relative risks to apply to each transition matrix in `x`.
#' @param index The indices of the transition probability matrices that `rr`  is applied to. 
#' Should either be a matrix where the first column denotes a transition probability matrix row 
#' and the second column denotes a transition probability matrix column or a list
#' where each element is a vector of length 2 with the first element denoting
#'  a transition probability matrix row and the second column denoting a transition
#'  probability matrix column.
#' @param complement Denotes indices of transition probability matrices that are
#' "complements" (i.e., should be computed as \eqn{1} less the sum of all other
#' elements in that row). Should be in the same format as `index`, If `NULL`, then
#' the diagonals of the matrix are assumed to be the complements. There can only be one
#' complement for each row in a transition probability matrix.
#' @export
apply_rr <- function(x, rr, index, complement = NULL) {
  
}