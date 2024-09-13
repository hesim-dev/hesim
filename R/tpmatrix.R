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
#' `data.table`, and `tpmatrix`. For `as_tbl2()` a 3-dimensional array 
#' where each slice is a square matrix.
#' @param output The class of the object returned by the function. Either 
#' a `data.table`, `data.frame`, `matrix`, or [`tpmatrix`].
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
#' as_tbl2(as_array3(pmat))
#' as_tbl2(as_array3(pmat), prefix = "p_", sep = ".")
#' @return For `as_array3()` a 3-dimensional array of square matrices; 
#' for `as_tbl2()` a 2-dimensional tabular object as specified by `output`.
#' 
#' @seealso [`tpmatrix`]
#' @export
as_array3 <- function(x) {
  if (length(dim(x)) != 2) stop("'x' must be a 2-dimensional object.")
  n_states <- sqrt(ncol(x))
  if (!is_whole_number(n_states)){
    stop("'x' must contain square matrices.", call. = FALSE)
  }
  y <- aperm(array(c(t(x)),
                   dim = c(n_states, n_states, nrow(x))),
             c(2, 1, 3))
  if (!is.null(attr(x, "states"))) {
    states <- attr(x, "states")
    dimnames(y) <- list(states, states, NULL)
  }
  return(y)
}

#' @name as_array3
#' @param prefix,sep Arguments passed to [tpmatrix_names()] for naming
#' the transition probability columns. The `states` argument is based on
#' the column names (i.e., names of the second dimension) of array; 
#' if `NULL`, then states are named `s1`, ..., `sh` where h is 
#' the number of states.
#' @export
as_tbl2 <- function(x, output = c("data.table", "data.frame", "matrix", "tpmatrix"),
                    prefix = "", sep = "_") {
  output <- match.arg(output)
  n_col <- dim(x)[1] * dim(x)[2]
  y <- matrix(c(aperm(x, c(2, 1, 3))), ncol = n_col, byrow = TRUE)
  if (is.null(colnames(x))) {
    states <- paste0("s", 1:dim(x)[1])
  } else{
    states <- colnames(x)
  }
  colnames(y) <- tpmatrix_names(states = states, prefix = prefix, sep = sep)
  if (output == "data.table") {
    return(as.data.table(y))
  } else if (output == "data.frame") {
    return(as.data.frame(y))
  } else if (output == "matrix") {
    return(y)
  } else {
    return(tpmatrix(y))
  }
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
#' `paste0(prefix, states[i], sep, states[j])`.
#' 
#' @param states A character vector of the names of health states in the 
#' transition matrix.
#' @param prefix A prefix that precedes the described transitions between states
#' used to name a transition. For example, if `prefix = "p_"` (and `sep = "_"`), 
#' then a transition between state `i` and state `j` will be of the form 
#' `"p_states[i]_states[j]"`; similarly, if `prefix = ""`, then the same 
#' transition will be named `"states[i]_states[j]"`. 
#' @param sep A character string to separate the terms representing
#' state `i` and state `j`. For instance, if `sep = "."`, the resulting name
#' will be of the form `"states[i].states[j]"`. 
#' @examples 
#' tpmatrix_names(LETTERS[1:4])
#' tpmatrix_names(LETTERS[1:4], prefix = "")
#' tpmatrix_names(LETTERS[1:4], prefix = "", sep = ".")
#' 
#' @return A character vector containing a name for each element of the transition
#' probability matrix encompassing all possible transitions. 
#' 
#' @seealso See `tpmatrix()`, which uses `tpmatrix_names()` to name the columns
#' of the returned object. 
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
  return(x[, start:(start + n_states - 1), drop = FALSE])
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
#' @param complement Either a character vector or a numeric vector denoting the 
#' transitions (i.e., the columns of the tabular object formed from `...`) that
#' are complementary (see "Details" below). If a character vector, each element
#' should be the name of a column in the tabular object; if a numeric vector,
#' each element should be the index of a column in the tabular object. 
#' 
#' @param states,prefix,sep Arguments passed to [tpmatrix_names()] for naming
#' the columns. If `states = NULL` (the default), then the states are named 
#' `s1`, ..., `sh` where `h` is the number of health states. 
#' 
#' @details A `tpmatrix` is a 2-dimensional tabular object that stores flattened
#' square transition probability matrices in each row. Each transition probability
#' matrix is filled rowwise. The complementary probability (equal to \eqn{1} 
#' minus the sum of the probabilities of all other elements in a row of a 
#' transition probability matrix) can be conveniently referred to as `C` or 
#' specified with the `complement` argument. There can only be one complement 
#' for each row in a transition probability matrix.
#' 
#' @return Returns a `tpmatrix` object that inherits from `data.table`
#' where each column is an element of the transition probability matrix with
#' elements ordered rowwise. 
#' 
#' 
#' @examples 
# Pass vectors
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
#' # Use the 'complement' argument
#' pmat <- data.frame(s1_s1 = 0, s1_s2 = .5, s2_s1 = .3, s2_s2 = 0)
#' tpmatrix(pmat, complement = c("s1_s1", "s2_s2"))
#' tpmatrix(pmat, complement = c(1, 4)) # Can also pass integers
#' 
#' # Can control column names
#' tpmatrix(pmat, complement = c(1, 4),
#'          states = c("state1", "state2"), sep = ".")
#' 
#' @seealso A `tpmatrix` is useful because it provides a convenient
#' way to construct a [`tparams_transprobs`] object, which is the object in
#' `hesim` used to specify the transition probabilities required to simulate 
#' Markov chains with the [`CohortDtstmTrans`] class. See the 
#' [`tparams_transprobs`] documentation for more details.
#' 
#' The [summary.tpmatrix()] method can be used to summarize a `tpmatrix` 
#' across parameter samples.
#' 
#' @export
tpmatrix <- function(..., complement = NULL, states = NULL,
                     prefix = "", sep = "_"){
  # Evaluate
  m_def <- define_tpmatrix(...)
  m <- eval_dots(m_def, as.list(parent.frame()))
  m_complement <- attr(m, "complement")
  m <- as.data.table(m)
  
  # Update complement
  if (is.null(complement)) {
    complement <- m_complement
  } else if (inherits(complement, c("numeric", "integer"))) {
    complement <- ifelse(1:ncol(m) %in% complement, 1, m_complement)
  } else if (inherits(complement, "character")) {
    complement <- ifelse(colnames(m) %in% complement, 1, m_complement)
  } else {
    stop("'complement' must either be a vector of integers or a character vector.")
  }

  # Some checks
  n_states <- sqrt(ncol(m))
  if (is.null(states)) states <- paste0("s", 1:n_states)
  if (!is_whole_number(n_states)){
    stop("tpmatrix() must be a square matrix.", call. = FALSE)
  }
  if (length(states) != n_states) {
    stop(paste0("The length of 'states' must equal the square root of the ",
                 "number of elements in the transition probability matrix."),
                call. = FALSE)
  }
  
  # Naming the columns
  colnames(m) <- tpmatrix_names(states = states, prefix = prefix, sep = sep)

  # Replace complement
  replace_C(m, complement)
  
  # Return
  setattr(m, "class", c("tpmatrix", "data.table", "data.frame"))
  setattr(m, "states", states)
  return(m)
}

#' Summarize transition probability matrix
#' 
#' Summarize a [`tpmatrix`] object storing transition probability matrices. 
#' Summary statistics are computed for each possible transition. 
#' 
#' @param object A [`tpmatrix`] object.
#' @param id A [`tpmatrix_id`] object for which columns contain the ID variables 
#' for each row in `object`. If not `NULL`, then transition probability matrices
#' are summarized by the ID variables in `id`. 
#' @param probs  A numeric vector of probabilities with values in `[0,1]` used
#' to compute quantiles. Computing quantiles can be slow when `object` is large,
#' so the default is `NULL`, meaning that no quantiles are computed.
#' @param unflatten If `FALSE`, then each column containing a summary statistic
#'  is a vector and the generated table contains one row 
#'  (for each set of ID variables) for each possible transition; if
#' `TRUE`, then each column stores a list of `matrix` objects containing
#' transition probability matrices formed by "unflattening" the one-dimensional
#' vectors. See "Value" below for additional details.
#' @param ... Additional arguments affecting the summary. Currently unused.
#' 
#' @return If `unflatten = "FALSE"` (the default), then a [`data.table::data.table`]
#' is returned with columns for (i) the health state that is being transitioned
#' from (`from`), (ii) the health state that is being transitioned to (`to`)
#' (iii) the mean of each parameter across parameter samples (`mean`),
#' (iv) the standard deviation of the parameter samples (`sd`), and
#' (v) quantiles of the parameter samples corresponding to the `probs` argument. 
#' 
#' If, on the other hand, `unflatten = "TRUE"`, then the parameters are unflattened
#' to form transition probability matrices; that is, the `mean`, `sd`, and 
#' quantile columns are (lists of) matrices. 
#' 
#' In both cases, if `id` is not `NULL`, then the ID variables are also 
#' returned as columns.
#' @examples 
#' library("data.table")
#' hesim_dat <-  hesim_data(strategies = data.table(strategy_id = 1:2),
#'                                 patients = data.table(patient_id = 1:3))
#' input_data <- expand(hesim_dat, by = c("strategies", "patients"))
#' 
#' # Summarize across all rows in "input_data"
#' p_12 <- ifelse(input_data$strategy_id == 1, .8, .6)
#' p <- tpmatrix(
#'   C, p_12,
#'   0, 1
#' )
#' 
#' ## Summary where each column is a vector
#' summary(p)
#' summary(p, probs = c(.025, .975))
#' 
#' ## Summary where each column is a matrix
#' ps <- summary(p, probs = .5, unflatten = TRUE)
#' ps
#' ps$mean
#' 
#' # Summarize by ID variables
#' tpmat_id <- tpmatrix_id(input_data, n_samples = 2) 
#' p_12 <- ifelse(tpmat_id$strategy_id == 1, .8, .6)
#' p <- tpmatrix(
#'   C, p_12,
#'   0, 1
#' )
#' 
#' ## Summary where each column is a vector
#' summary(p, id = tpmat_id)
#'
#' ## Summary where each column is a matrix
#' ps <- summary(p, id = tpmat_id, unflatten = TRUE)
#' ps
#' ps$mean
#' @export
summary.tpmatrix <- function(object, id = NULL, probs = NULL, 
                             unflatten = FALSE, ...) {
  
  states <- attr(object, "states")
  n_states <- length(states)
  n_trans <- ncol(object)
  n_id <- nrow(object)
  
  # Convert transition probabilities to long format
  object_dt <- data.table(
    from = rep(rep(states, each = n_states), n_id),
    to = rep(rep(states, times = n_states), n_id),
    prob = c(t(object))
  )
  
  # Add ID variables if they are passed by user
  if (!is.null(id)) {
    check_is_class(id, "tpmatrix_id", "tpmatrix_id")
    
    # Get relevant ID variables
    id_cols <- attr(id, "id_vars")
    if ("time_start" %in% colnames(id)) id_cols <- c(id_cols, "time_start") 
    if ("time_stop" %in% colnames(id)) id_cols <- c(id_cols, "time_stop") 
    
    # Lengthen ID table to same length as long transition probabilities
    # and column bind relevant ID variables
    id2 <- id[rep(seq_len(nrow(id)), each = n_trans)]
    id2 <- id2[, id_cols, with = FALSE]
    object_dt <- cbind(id2, object_dt)
  } else {
    by <- NULL
  }
  
  # Summarize and return
  summarize_transprobs_dt(object_dt, probs = probs, unflatten = unflatten,
                          states = states)
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
#' @return Returns a `tpmatrix_id` object that inherits from `data.table` with 
#' the same columns in `object` repeated `n_samples` times. That is, to facilitate
#' creation of a `tparams_transprobs()` object,  there is one row for each 
#' parameter sample, treatment strategy, patient, and optionally time interval.
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
    stop(paste0("'object' must be expanded by 'strategy_id', 'patient_id', ", 
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
#' A generic function for creating transition intensity matrices where 
#' elements represent the instantaneous risk of moving between health states.
#' @param x An `R` object.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @keywords internal
#' @export
qmatrix <- function (x, ...) {
  UseMethod("qmatrix", x)
}

#' Transition intensity matrix from tabular object
#' 
#' Creates transition intensity matrices where elements represent
#' the instantaneous risk of moving between health states. 
#' 
#' @param x A two-dimensional tabular object containing
#' elements of the transition intensity matrix. A column represents a transition
#' from state \eqn{r} to state \eqn{s}. Each row represents elements of a different
#' transition intensity matrix. See "Details" for more information.
#' @param trans_mat Just as in [`IndivCtstmTrans`], a transition matrix 
#' describing the states and transitions in a multi-state model.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' 
#' @details The object `x` must only contain non-zero and non-diagonal elements
#' of a transition intensity matrix. The diagonal elements are automatically computed
#' as the negative sum of the other rows. 
#' 
#' @return An array of transition intensity matrices with the third dimension 
#' equal to the number of rows in `x`.
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
#' @seealso [qmatrix.msm()]
#' @export
qmatrix.matrix <- function(x, trans_mat, ...){
  q <- as.matrix(x)
  trans <- c(t(trans_mat))
  n_states <- nrow(trans_mat)
  qmat <- matrix(0, nrow = nrow(q), ncol = length(trans))
  qmat[, which(!is.na(trans))] <- q
  colnames(qmat) <- tpmatrix_names(states = paste0("s", 1:n_states),
                                   prefix = "")
  qmat <- replace_Qdiag(qmat, n_states)
  return(as_array3(qmat))
}

#' @rdname qmatrix.matrix
#' @export
qmatrix.data.table <- function(x, trans_mat, ...) {
  return(qmatrix(as.matrix(x), trans_mat))
}

#' @rdname qmatrix.matrix
#' @export
qmatrix.data.frame <- function(x, trans_mat, ...) {
  return(qmatrix(as.matrix(x), trans_mat))
}

#' Transition intensity matrix from `msm` object
#' 
#' Draw transition intensity matrices for a probabilistic sensitivity analysis
#' from a fitted `msm` object.
#' 
#' @param x A [`msm::msm`] object.
#' @param newdata A data frame to look for variables with which to predict. A 
#' separate transition intensity matrix is predicted based on each row in
#' `newdata`. Can be `NULL` if no covariates are included in the fitted `msm`
#' object.
#' @param uncertainty Method used to draw transition intensity matrices. If `"none`",
#' then point estimates are used. If `"normal"`, then samples are drawn from the 
#' multivariate normal distribution of the regression coefficients. 
#' @param n Number of random observations of the parameters to draw.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' 
#' @return An array of transition intensity matrices with the third dimension 
#' equal to the number of rows in `newdata`.
#' 
#' @examples 
#' library("msm")
#' set.seed(101)
#'  qinit <- rbind(
#'    c(0, 0.28163, 0.01239),
#'    c(0, 0, 0.10204),
#'    c(0, 0, 0)
#'  )
#' fit <- msm(state_id ~ time, subject = patient_id, 
#'            data = onc3p[patient_id %in% sample(patient_id, 100)],
#'            covariates = list("1-2" =~ age + strategy_name), 
#'            qmatrix = qinit)
#' qmatrix(fit, newdata = data.frame(age = 55, strategy_name = "New 1"),
#'         uncertainty = "none")
#' qmatrix(fit, newdata = data.frame(age = 55, strategy_name = "New 1"),
#'         uncertainty = "normal",  n = 3)
#' 
#' @seealso `qmatrix.matrix()`
#' @export
qmatrix.msm <- function(x, newdata = NULL, uncertainty = c("normal", "none"), n = 1000, 
                        ...) {
  uncertainty <- match.arg(uncertainty)
  if (is.null(newdata) && x$qcmodel$ncovs > 0) {
    stop("'newdata' cannot be NULL if covariates are included in 'x'.")
  }
  which_base <- which(x$paramdata$plabs=="qbase")
  which_cov <- which(x$paramdata$plabs=="qcov") 

  # Simulate distribution of parameters
  if (uncertainty == "normal") {
    beta <- normboot(x, n)
  } else if (uncertainty == "none"){
    beta <- matrix(x$paramdata$params[c(which_base, which_cov)], nrow = 1)
  }
  
  # Create model matrix from newdata
  if (is.null(newdata)) {
    X <- matrix(1, nrow = 1, ncol = 1)
  } else{
    mfo <- stats::model.frame(x$covariates, x$data$mf)
    tt <- attr(mfo, "terms")
    Terms <- stats::delete.response(tt)
    mf <- stats::model.frame(Terms, newdata, xlev = stats::.getXlevels(tt, mfo))
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      stats::.checkMFClasses(cl, mf)
    X <- stats::model.matrix(Terms, mf)
    if (x$center) {
      X <- cbind(X[, 1, drop = FALSE],
                 sweep(X[, -1, drop = FALSE], 2, x$qcmodel$covmeans))
    } 
  }
  
  # Predict each element of the transition intensity matrix
  res <- array(NA, dim = c(x$qmodel$nstates, x$qmodel$nstates, nrow(X) * nrow(beta)))
  tmat <- x$qmodel$imatrix
  tmat[tmat == 0] <- NA
  n_params <- sum(!is.na(tmat))
  n_covs <- x$qcmodel$ncovs
  q <- matrix(0, nrow = dim(res)[3], ncol = n_params)
  for (i in 1:n_params) {
    ind <- seq(from = i, by = n_params, length.out = n_covs + 1)
    q[, i] <- c(X %*% t(beta[, ind, drop = FALSE]))
  }
  return(qmatrix(exp(q), tmat))
}

# Matrix exponential -----------------------------------------------------------
#' Matrix exponential
#' 
#' This is a wrapper around [`msm::MatrixExp()`] that computes the exponential
#' of multiple square matrices. 
#' 
#' @param x An array of matrices. 
#' @param t An optional scaling factor for `x`. 
#' @param ... Arguments to pass to [`msm::MatrixExp()`].
#' 
#' @details This function is most useful when exponentiating transition intensity
#' matrices to produce transition probability matrices. To create transition
#' probability matrices for discrete time state transition models with annual
#' cycles, set `t=1`. An array of matrices is returned which can be used
#' to create the `value` element of a [`tparams_transprobs`] object. See 
#' [`qmatrix()`] for an example. 
#' 
#' @return Returns an array of exponentiated matrices. If `length(t) > 1`, then
#' `length(t)` arrays are returned for each element in `x`.
#' 
#' @seealso [`qmatrix.msm()`], [`qmatrix.data.table()`]
#' @export
expmat <- function(x, t = 1, ...) {
  if (!is.array(x)) {
    stop("'x' must be an array.")
  }
  if (length(dim(x)) == 2) {
    x <- array(x, dim = c(dim(x)[1], dim(x)[1], 1))
  }
  n <- dim(x)[3]
  n_states <- dim(x)[1]
  n_t <- length(t)
  res <- array(NA, dim = c(n_states, n_states, n * n_t))
  counter <- 1
  for (i in 1:n) {
    for (j in seq_along(t)) {
      res[,, counter] <- msm::MatrixExp(x[,, i], t = t[j], ...)
      counter <- counter + 1
    }
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
#' a transition probability matrix row and the second column denoting a transition
#' probability matrix column.
#' @param complement Denotes indices of transition probability matrices that are
#' "complements" (i.e.,  computed as \eqn{1} less the sum of all other
#' elements in that row). Should be in the same format as `index`. There can be
#' at most one complementary column in each row of a transition probability 
#' matrix. If `NULL`, then the diagonals are assumed to be the complements.
#' 
#' @details This function is useful for applying relative treatment effects measured
#' using relative risks to an existing transition probability matrix. For example,
#' a transition probability matrix for the reference treatment strategy may exist or
#' have been estimated from the data. Relative risks estimated from a meta-analysis
#' or network meta-analysis can then be applied to the reference transition probability
#' matrix. If the number of rows in `rr` exceeds `x`, then the arrays in `x` are 
#' recycled to the number of rows in `rr`, which facilitates the application of 
#' relative risks from multiple treatment strategies to a reference treatment.
#' 
#' @return A 3-dimensional array where each slice contains matrices of the same 
#' dimension as each matrix in `x` and the number of slices is equal to the number
#' of rows in `rr`.
#'  
#' @examples
#' p_12 <- c(.7, .5)
#' p_23 <- c(.1, .2)
#' x <- as_array3(tpmatrix(
#'   C, p_12, .1,
#'   0, C,     p_23,
#'   0, 0,     1
#' ))
#' 
#' # There are the same number of relative risk rows and transition probability matrices
#' rr_12 <- runif(2, .8, 1)
#' rr_13 <- runif(2, .9, 1)
#' rr <- cbind(rr_12, rr_13)
#' apply_rr(x, rr, 
#'          index = list(c(1, 2), c(1, 3)),
#'          complement = list(c(1, 1), c(2, 2)))
#'          
#' # There are more relative risk rows than transition probability matrices
#' rr_12 <- runif(4, .8, 1)
#' rr_13 <- runif(4, .9, 1)
#' rr <- cbind(rr_12, rr_13)
#' apply_rr(x, rr, 
#'          index = list(c(1, 2), c(1, 3)),
#'          complement = list(c(1, 1), c(2, 2)))
#' @export
apply_rr <- function(x, rr, index, complement = NULL) {
  n_rows <- dim(x)[1]
  if (is.list(index)) index <- do.call("rbind", index)
  if (is.null(complement)) {
    complement <- cbind(1:n_rows, 1:n_rows)
  } else {
    if (is.list(complement)) complement <- do.call("rbind", complement)
  }
  
  # Checks
  if (nrow(index) != ncol(rr)) {
    stop(paste0("'index' must contain the same number of matrix elements as the ",
                "number of columns in 'rr'."), call. = FALSE)
  }
  if (nrow(complement) > n_rows) {
    stop(paste0("The number of matrix elements in 'complement' cannot be larger than the ",
                "number of rows in 'x'."), call. = FALSE)
  }
  if (any(table(complement[, 1]) != 1)) {
    stop("There can only be one complementary column in each row.", call. = FALSE)
  }
  
  # Run C function
  C_apply_rr(x, rr, index - 1, complement - 1)
}