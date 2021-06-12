# Random number generation -----------------------------------------------------
#' Define and evaluate random number generation expressions
#' 
#' Random number generation expressions are used to 
#' randomly sample model parameters from suitable distributions for probabilistic
#' sensitivity analysis. These functions are typically used when evaluating
#' an object of class `model_def` defined using [define_model()].
#' 
#' @param expr An expression used to randomly draw variates for each parameter of
#' interest in the model. [Braces][base::Paren] should be used so that the result
#' of the last expression within the braces is evaluated. The expression must
#'  return a list where each element is either a `vector` or tabular object
#'  ( `matrix`, `data.frame`, or `data.table`). The length of the vector must
#'  either be or `n` and the number of rows in the tabular object must be `n`.
#' @param n Number of samples of the parameters to draw.
#' @param ... Additional arguments to pass to the environment used to evaluate
#' `expr`.
#' 
#' @details `hesim` contains a number of random number generation functions
#' that return parameter samples in convenient formats
#' and do not typically require the number of samples, `n`, as arguments 
#' (see [`rng_distributions`]). The random number generation expressions
#' are evaluated using `eval_rng()` and used within `expr`
#' in `define_rng()`. If a multivariate object is returned by `eval_rng()`,
#' then the rows are random samples and columns are 
#' distinct parameters (e.g., costs for each health state, elements of a 
#' transition probability matrix).  
#' 
#' @return `define_rng()` returns an object of class `rng_def`,
#' which is a list containing the unevaluated random number generation
#' expressions passed  to `expr`, `n`, and any additional arguments passed to 
#' `...` . `eval_rng()` evaluates the `rng_def` object and 
#' returns an `eval_rng` object containing the evaluated expression.
#' 
#' @seealso Parameters can be conveniently sampled from probability distributions 
#' using a number of random number generation functions (see [`rng_distributions`]). 
#' An economic model can be created with [create_CohortDtstm()] by using 
#' `define_rng()` (or a previously evaluated `eval_rng` object)  
#' alongside [define_tparams()] to define a model with [define_model()].
#' It can be useful to summarize an evaluated expression with `summary.eval_rng()`.
#' @examples  
#' params <- list(
#'   alpha = matrix(c(75, 25, 33, 67), byrow = TRUE, ncol = 2),
#'   inptcost_mean = c(A = 900, B = 1500, C = 2000),
#'   outptcost_mean = matrix(c(300, 600, 800,
#'                             400, 700, 700),
#'                            ncol = 3, byrow = TRUE)
#' )
#' rng_def <- define_rng({
#'   aecost_mean <- c(500, 800, 1000) # Local object not 
#'                                    # not returned by eval_rng()
#'   list( # Sampled values of parameters returned by eval_rng()
#'     p = dirichlet_rng(alpha), # Default column names
#'     inptcost = gamma_rng(mean = inptcost_mean, # Column names based on 
#'                          sd = inptcost_mean),  # named vector
#'     outptcost = outptcost_mean, # No column names because
#'                                 # outptcost_mean has none.
#'     aecost = gamma_rng(mean = aecost_mean, # Explicit naming of columns
#'                        sd = aecost_mean,
#'                        names = aecost_colnames)
#'   )
#' }, n = 2, aecost_colnames = c("A", "B", "C")) # Add aecost_colnames to environment
#' params_sample <- eval_rng(x = rng_def, params)
#' summary(params_sample)
#' params_sample
#' @export
define_rng <- function(expr, n = 1, ...){
  x <- c(list(expr = substitute(expr), n = n), list(...))
  class(x) <- "rng_def"
  return(x)
}

#' @param x An object of class `rng_def` created with `define_rng()`.
#' @param params A list containing the values of parameters for random number 
#' generation. Each element of the list should either be a `vector`,
#'  `matrix`, `data.frame`, or `data.table`
#' @param check Whether to check the returned output so that (i) it returns a list
#' and (ii) each element has the correct length or number of rows. Default is `FALSE`,
#' meaning that any output can be returned. This is always `TRUE` when used inside
#' [define_model()]. 
#' @export
#' @rdname define_rng
eval_rng <- function(x, params = NULL, check = FALSE){
  y <- eval(x$expr, envir = c(x[-1], params)) # -1 is the position of "expr"
  if (!inherits(y, "list")){
    stop("define_rng() must return a list.", call. = FALSE)
  } 
  class(y) <- "eval_rng"
  if (check) check(y)
  attr(y, "n") <- x$n
  attr(y, "checked") <- check
  return(y)
}

#' @export
as.list.eval_rng <- function(x, ...) {
  class(x) <- "list"
  x
}

#' @export
c.eval_rng <- function(...) {
  dots <- list(...)
  n <- attr(dots[[1]], "n")
  class(dots[[1]]) <- "list"
  x <- do.call("c", dots)
  class(x) <- "eval_rng"
  attr(x, "n") <- n
  x
}

#' Sumnmarize `eval_rng` object
#' 
#' Summarize the model parameters randomly sampled for probabilistic sensitivity 
#' analysis with [eval_rng()]. 
#' @param object,x An [`eval_rng`] object.
#' @param probs A numeric vector of probabilities with values in `[0,1]` used to 
#' compute quantiles with [stats::quantile()].
#' @param sep When a list element returned by `eval_rng` is a tabular object,
#' the parameter name is created by concatenating the name of the list element 
#' with the columns of the tabular object. The `sep` argument determines the 
#' character string used to separate the terms.
#' @param ... For the print method, arguments to pass to `summary.eval_rng()`. 
#' 
#' @return `summary.eval_rng()` returns a [`data.table`] with columns for
#' (i) the name of the parameter (`param`), (ii) the mean of the parameter
#' samples (`mean`), (iii) the standard deviation of the parameter samples (`sd`),
#' and (iv) quantiles of the parameter samples corresponding
#' to the `probs` argument. `print.eval_rng()` prints the output of 
#' `summary.eval_rng()` to the console. 
#' 
#' @seealso See [eval_rng()] for an example. 
#' @export
summary.eval_rng <- function(object, probs = c(.025, .975), sep = "_",  ...) {
  
  apply_quantile <- function(x, probs) {
    y <- apply(x, 2, stats::quantile, probs = probs)
    if (length(probs) == 1) {
      y <- matrix(y, ncol = 1)
      colnames(y) <- paste0(probs * 100, "%")
      return(y)
    } else{
      return(t(y))
    }
  }
  
  fun <- function(x, name, sep) {
    if (is_1d_vector(x)) {
      as.data.table(t(c(
        param = name,
        mean = mean(x),
        sd = stats::sd(x),
        stats::quantile(x, probs = probs)
      )))
    } else {
      p <- if (!is.null(colnames(x))) colnames(x) else paste0("v", 1:ncol(x))
      data.table(
        param = paste0(name, sep, p),
        mean = apply(x, 2, mean),
        sd = apply(x, 2, stats::sd),
        apply_quantile(x, probs)
      )
    }
  }
  
  res_list <- lapply(seq_along(object), function(i, x, names) {
    fun(x[[i]], names[i], sep)
  }, x = object, names = names(object))
  rbindlist(res_list)
}

#' @rdname summary.eval_rng
#' @export
print.eval_rng <- function(x, ...) {
  cat("A summary of the \"eval_rng\" object:")
  cat("\n\n")
  print(summary(x, ...))
  invisible(x)
}

check.eval_rng <- function(object){
  
  object <- as.list(object)
  
  # Number of samples
  fun <- function(z){
    if (length(dim(z)) == 2){
      return(list(n = nrow(z), is_vector = FALSE))
    } else if (is_1d_vector(z)) {
      return(list(n = length(z), is_vector = TRUE))
    } else {
      stop(paste0("Each element returned by define_rng() must either be a ",
                  "vector or a tabular object."))
    }
  }  
  
  n_elem <- length(object)
  is_vector <- n <- rep(NA, n_elem)
  for (i in 1:n_elem){
    res <- fun(object[[i]])
    is_vector[i] <- res$is_vector
    n[i] <- res$n
  }
  if(any(!(n == n[1] | (is_vector == TRUE & n == 1)))){
    stop(paste0("The number of samples produced by define_rng() must be ",
                "equal to n unless a scalar (of length 1) is returned."),
         call. = FALSE)
  }
}

rng_colnames <- function(old_names, params, new_names = NULL){
  k <- length(params[[1]])
  if (is.null(new_names)){ # Case where names == NULL
    # Get names of first non NULL named parameter vector
    vec_names <- NULL
    for (i in 1:length(params)){
      vec_names_i <- names(params[[i]])
      if (!is.null(vec_names_i)){
        vec_names <- vec_names_i
        break
      }
    }
    if (!is.null(vec_names)){
      new_names <- vec_names
    } else{
      new_names<- paste0("v", as.character(1:k))
    }
  } 
  return(new_names)
}

mom_fun_rng <- function(n, rng_fun, mom_fun, mean, sd){
  which_fixed <- which(sd == 0)
  which_rng <- which(sd != 0)
  if (length(which_rng) > 0){
    rng_params <- do.call(mom_fun, list(mean[which_rng], sd[which_rng]))
    x <- do.call(rng_fun, c(list(n = n * length(which_rng)),
                            rng_params))
  }
  if (length(which_fixed) > 0){ # Use fixed() when sd == 0
    x_fixed <- rep(mean[which_fixed], length = n * length(which_fixed))
    if (length(which_rng) > 0){ # Combine fixed() and rng()
      x <- cbind(matrix(x, byrow = TRUE, nrow = n),
                 matrix(x_fixed, byrow = TRUE, nrow = n))
      x <- x[, order(c(which_rng, which_fixed))]
      x <- c(t(x))
    } else{ # Only used fixed()
      x <- x_fixed
    }
  }
  return(x)
}

#' Generate variates for univariate distributions
#' @keywords internal
uv_rng <- function(n, params, rng_fun, mom_fun = NULL, names = NULL){
  
  # Check length of all parameter vectors are equal
  k <- sapply(params, length)
  if (!all(k == k[1])){
    stop("The length of all vectors of parameters (e.g., mean, sd) must be equal.")
  }
  k <- k[1]
  
  # Random number generation
  if (!is.null(mom_fun)){
    x <- mom_fun_rng(n, rng_fun, mom_fun, params[[1]], params[[2]])
  } else{
    x <- do.call(rng_fun, c(list(n = n * k), params))
  }
  
  ## If k > 1 convert to matrix
  if (k > 1){
    x <- matrix(x, byrow = TRUE, nrow = n)
  } 
  
  # Make nice names if a matrix
  if (is.matrix(x)){
    colnames(x) <- rng_colnames(colnames(x), params, names)
  }
  
  # Return
  if (k == 1){
    return(x) # Return vector
  } else{
    return(data.table(x)) # Return data.table
  }
}

#' Random number generation distributions
#' 
#' A collection of functions for randomly generating deviates from probability
#' distributions with [define_rng()]. 
#' @param mean,sd Mean and standard deviation of the random variable.
#' @param names Names for columns if an object with multiple columns is returned 
#' by the function. 
#' @param n The number of random samples of the parameters to draw. Default is
#' the value of `n` in the environment in which the function is called, which 
#' can be useful when used inside `define_rng` because it means that a value does
#' not need to be explictly passed to `n`.
#' @param ... Additional arguments to pass to underlying random number generation
#' functions. See "details". 
#' 
#' @details These functions are not exported and are meant for use with
#' [define_rng()]. They consequently assume that the number of samples to draw, `n`,
#' is defined in the parent environment. Convenience random number generation 
#' functions include:
#' \describe{
#' \item{`beta_rng()`}{If `mean` and `sd` are both not `NULL`, then
#'  parameters of the beta distribution are derived using
#' the methods of moments with [mom_beta()]. Beta variates are generated with
#'  [stats::rbeta()].}
#' \item{`custom()`}{Use previously sampled values from a custom probability distribution.
#' There are three possibilities: (i) if `n` is equal to the number previously 
#' sampled values (say `n_samples`), then `x` is returned as is; (ii) if
#'  `n` < `n_samples`, then samples from `x` are sampled without replacement;
#' and (iii) if `n` > `n_samples`, then samples from `x` are sampled with replacement
#' and a warning is provided.}
#' \item{`dirichlet_rng()`}{Dirichlet variates for each row in the matrix are
#' generated with [rdirichlet_mat()]. The sampled values are stored in a `data.table`
#'  where there is a column for each element of `alpha` 
#'  (with elements ordered rowwise).}
#' \item{`fixed()`}{This function should be used when values of the variable
#'  of interest are fixed (i.e., they are known with certainty). If `length(est) > 1`, 
#'  an `n` by `length(est)` `data.table` is returned meaning that each element of `est`
#'  is repeated `n` times; otherwise (if `length(est) == 1`), a vector is returned
#'  where `est` is repeated `n` times.}
#' \item{`gamma_rng()`}{The parameters of the gamma distribution are derived using
#'  the methods of moments with [mom_gamma()] and gamma variates are generated 
#'  with [stats::rgamma()].}
#' \item{`lognormal_rng()`}{Lognormal variates are generated with [stats::rlnorm()].}
#' \item{`multi_normal_rng()`}{Multivariate normal variates are generated with [MASS::mvrnorm()].}
#' \item{`normal_rng()`}{Normal variates are generated with [stats::rnorm()].}
#' \item{`uniform_rng()`}{Uniform variates are generated with [stats::runif()].}
#' 
#' }
#' 
#' @return
#' Functions either return a vector of length `n` or an `n` by `k` `data.table`.
#'  Multivariate distributions always return a `data.table`. If a 
#' univariate distribution is used, then a `data.table` is returned if each 
#' parameter is specified as a vector with length greater than 1; otherwise, if
#'parameters are scalars, then a vector is returned. In the `data.table` case,
#'  `k` is equal to the length of the parameter vectors 
#' entered as  arguments. For example, if the probability distribution contained
#'  `mean` as an argument and `mean` were 
#' of length 3, then an `n` by 3 matrix would be returned. The length of all 
#' parameter vectors must be the same. For instance, if the vector `mean` 
#' were of length 3 then all additional parameters (e.g., `sd`)
#' must also be of length 3.
#' 
#' If a `data.table` is returned by a distribution, then its column names are set
#'  according to the following hierarchy:
#' \enumerate{
#'   \item With the `names` argument if it is not `NULL`
#'   \item With the names of the parameter vectors if they are named vectors. If there
#'   are multiple parameter vector arguments, then the names of the first parameter
#'   vector with non `NULL` names is used. For instance, if `mean` and `sd` are
#'   both arguments to a random number generation function and `mean` is a
#'    named vector, then the names from the vector `mean` are used.
#'   \item As `v1`, ..., `vk` if the `names` argument is `NULL` and there are no
#'   named parameter vectors.
#' }
#' 
#' 
#' @name rng_distributions
#' @seealso [define_rng()]
NULL

#' @param shape1,shape2 Non-negative parameters of the Beta distribution.
#' @name rng_distributions
beta_rng <- function(shape1 = 1, shape2 = 1,
                     mean = NULL, sd = NULL, names = NULL,
                     n = parent.frame()$n){
  if (!is.null(mean) & !is.null(sd)){
    return(uv_rng(n = n, 
                  params = list(mean, sd) ,
                  rng_fun = "rbeta",
                  mom_fun = "mom_beta",
                  names = names))  
  } else{
    return(uv_rng(n = n, 
                  params = list(shape1, shape2) ,
                  rng_fun = "rbeta",
                  names = names)) 
  }
}

#' @param alpha A matrix where each row is a separate vector of shape parameters.
#' @rdname rng_distributions
dirichlet_rng <- function(alpha, names = NULL, n = parent.frame()$n){
  x <- rdirichlet_mat(n = n, alpha = alpha, output = "data.table")
  
  # Column names
  make_names <- function(x, y){
    y <- paste0("_", y)
    return(c(t(outer(x, y, paste0))))
  }
  if (is.null(names)){
    if (!is.null(rownames(alpha)) & !is.null(colnames(alpha))){
      colnames(x) <- make_names(rownames(alpha), colnames(alpha))
    } else{
      colnames(x) <- make_names(paste0("s", 1:nrow(alpha)), 
                                paste0("s", 1:ncol(alpha)))
    }
  } else{
    colnames(x) <- names
  }
  
  # Return
  return(x)
}

#' @param est A vector of estimates of the variable of interest. 
#' @rdname rng_distributions
fixed <- function(est, names = NULL, n = parent.frame()$n){
  stopifnot(is.vector(est))
  if (length(est) == 1){
    return(rep(est, n))
  } else{
    x <- matrix(est, byrow = TRUE, nrow = n,
                ncol = length(est))
    colnames(x) <- rng_colnames(colnames(x), list(est), names)
    return(data.table(x))
  }
}

#' @param x A numeric `vector`, `matrix`, `data.frame`, or `data.table` containing 
#' random samples of the variable of interest from a suitable probability distribution. This would 
#' typically be a posterior distribution from a Bayesian analysis. 
#' @rdname rng_distributions
custom <- function(x, names = NULL, n = parent.frame()$n){
  stopifnot(is.numeric(x))
  n_dims <- length(dim(x)) 
  if (n_dims > 2){
    stop("'x' must either be a vector, matrix, data.frame. or data.table.")
  }
  if (n_dims < 2){
    x <- matrix(x, ncol = 1)
  }

  # Return samples from posterior distribution
  samples <- sample_from_posterior(n = n,
                                   n_samples = nrow(x))
  x <- x[samples, ]
  
  # Return
  if (n_dims < 2){
    return(c(x))
  } else{
    if (!is.null(names)){
      colnames(x) <- names
    } else if (is.null(colnames(x))){
      colnames(x) <- paste0("v", as.character(1:ncol(x)))
    }
    return(data.table(x))
  }
}

#' @rdname rng_distributions
gamma_rng <- function(mean, sd, names = NULL, n = parent.frame()$n){
  return(uv_rng(n = n, 
                params = list(mean, sd),
                rng_fun = "rgamma",
                mom_fun = "mom_gamma",
                names = names))
}

#' @param meanlog,sdlog Mean and standard deviation of the distribution on the 
#' log scale.
#' @rdname rng_distributions
lognormal_rng <- function(meanlog, sdlog, names = NULL, n = parent.frame()$n){
  return(uv_rng(n = n, 
                params = list(meanlog, sdlog),
                rng_fun = "rlnorm",
                names = names))
}

#' @param mu,Sigma `mu` is a vector giving the means of the variables and
#' `Sigma` is a positive-definite symmetric matrix specifying the 
#' covariance matrix of the variables.
#' @rdname rng_distributions
multi_normal_rng <- function(mu, Sigma, names = NULL, n = parent.frame()$n, ...){
  m <- MASS::mvrnorm(n, mu = mu, Sigma = Sigma, ...)
  if (n == 1) {
    if (length(m) == 1) { # Case where n = 1 and a scalar is returned
      return(m)
    } else { # Otherwise convert to matrix
      m <- matrix(m, nrow = 1)
    }
  }
  
  colnames(m) <- rng_colnames(colnames(m), list(mu, Sigma), names)
  return(data.table(m)) 
}

#' @rdname rng_distributions
normal_rng <- function(mean, sd, names = NULL, n = parent.frame()$n){
  return(uv_rng(n = n, 
                params = list(mean, sd),
                rng_fun = "rnorm",
                names = names))
}

#' @param min,max Lower and upper limits of the distribution. Must be finite.
#' @rdname rng_distributions
uniform_rng <- function(min, max, names = NULL, n = parent.frame()$n){
  return(uv_rng(n = n, 
                params = list(min, max),
                rng_fun = "runif",
                names = names))
}

# Transformed parameters -------------------------------------------------------
#' Define and evaluate transformed parameter expressions
#' 
#' Transformed parameter expressions are used to transform the parameter
#' values sampled with [eval_rng()] as a function of input data 
#' (treatment strategies and patients) and time
#' intervals. These functions are used when evaluating an object of class
#' `model_def` defined using [define_model()]. The transformed parameters 
#' are ultimately converted into [tparams] objects and used to simulate outcomes with an
#' economic model.
#' @param expr Expressions used to transform parameters. As with [define_rng()], 
#' [braces][base::Paren] should be used so that the result
#' of the last expression within the braces is evaluated. The expression 
#' must return a named list with the following possible elements:
#' * *tpmatrix*: The transition probability matrix used to simulate transition
#' probabilities in the economic model. This should either be the output of
#' [tpmatrix()] or a 3-dimensional array as in [tparams_transprobs()].
#' * *utility*: The utility values to attach to states and used to simulate
#' quality-adjusted life-years in the economic model. Either a vector (in
#' which case utility is the same in each health state) or a
#'  `data.table`/`data.frame`/`matrix` with a column for each (non-death)
#'  health state.
#'  * *costs*:  A named list of costs for each category used to simulate
#'  costs in the economic model. Each element of the
#'  list must be in the same format as `utility`.
#'     
#' @param times Distinct times denoting the stopping time of time intervals.
#' @param ... Additional arguments to pass to the environment used to evaluate
#' `expr`.
#' 
#' @details `define_tparams()` is evaluated when creating economic models as a 
#' function of `model_def` objects defined with [define_model()]. Operations 
#' are "vectorized" in the sense that they are performed for each unique combination
#' of `input_data` and `params`. `expr` is evaluated in an environment including
#' each variable from `input_data`, all elements of `rng_params`, and a variable
#' `time` containing the values from `times`. The `time` variable can be used
#' to create models where parameters vary as a function of time. 
#' `eval_tparams()` is not exported and is only meant for use within [eval_model()].
#' 
#' @return [define_tparams()] returns an object of class `tparams_def`,
#'  which is a list containing the unevaluated "transformation" expressions
#'  passed  to `expr`, `times`, and any additional arguments passed to 
#'   `...` . [eval_tparams()] evaluates the `tparams_def` object 
#'  and should return a list of transformed parameter objects.
#' @seealso [define_model()], [define_rng()]
#' @export
define_tparams <- function(expr, times = NULL, ...){
 if (is.null(times)) times <- 0
  x <- c(list(expr = substitute(expr), times = times), list(...))
  class(x) <- "tparams_def"
  return(x)
}

check_eval_tparams <- function(object){
  # Must be a list
  if (!inherits(object, "list")){
    stop("define_tparams() must return a list", call. = FALSE)
  } 
  if (!inherits(object[["costs"]], "list") & !is.null(object[["costs"]])){
    stop("The 'costs' element returned by define_tparams() must be a list", call. = FALSE)
  } 
  
}

#' @param x An object of class `tparams_def`.
#' @param input_data An object of class [expanded_hesim_data][expand.hesim_data()] (as
#' in [eval_model()]) expanded by the distinct times in `times`. 
#' @param rng_params Random samples of the parameters returned by [eval_rng()]. 
#' @rdname define_tparams
eval_tparams <- function(x, input_data, rng_params){
  x <- eval(x$expr, envir = c(x[-1], as.list(input_data), rng_params)) # -1 is the position of "expr"
  check_eval_tparams(x)
  return(x)
}

# Model definition -------------------------------------------------------------
#' Define and evaluate model expression
#' 
#' A model expression is defined by specifying random number generation functions 
#' for a probabilistic sensitivity analysis (PSA) and transformations of the sampled 
#' parameters as a function of `input_data`. The unevaluated expressions
#' are evaluated with `eval_model()` and used to generate the model inputs needed to 
#' create an economic model.
#' 
#' @param tparams_def A [tparams_def][define_tparams()] object or a list of 
#' [tparams_def][define_tparams()] objects. A list might be considered if time intervals
#' specified with the `times` argument in [define_tparams()] vary across parameters. 
#' Parameters for a transition probability matrix (`tpmatrix`), utilities (`utility`),
#' and/or cost categories (`costs`) are returned as a named list (see [define_tparams()]
#' for more details).
#' @param rng_def A [rng_def][define_rng()] object used to randomly draw samples
#' of the parameters from suitable probability distributions.
#' @param params 	Either (i) a list containing the values of parameters for random 
#' number generation or (ii) parameter samples that have already been randomly
#' generated using `eval_rng()`. In case (ii), `rng_def` should be `NULL`. 
#' @param n_states The number of health states (inclusive of all health states
#' including the the death state) in the model. If `tpmatrix` is 
#' an element returned by `tparams_def`, then it will be equal to the number of states 
#' in the transition probability matrix; otherwise it must be specified as an argument.
#'  
#' @details `eval_model()` evaluates the expressions in an object of class
#'  `model_def` returned by `define_model()` and is, in turn, used within 
#'  functions that instantiate economic models (e.g., [create_CohortDtstm()]).
#'  The direct output of `eval_model()` can also be useful for understanding and debugging 
#'  model definitions, but it is not used directly for simulation.
#'  
#'  Economic models are constructed as a function of input data and parameters:
#'  1. *Input data*: Objects of class [expanded_hesim_data][expand.hesim_data()] 
#'  consisting of the treatment strategies and patient population.
#'  2. *Parameters*: The underlying parameter estimates from the literature 
#'  are first stored in a list (`params` argument). Random number generation
#'  is then used to sample the parameters from suitable probability distributions
#'  for the PSA (`rng_def` argument). Finally, the sampled parameters are
#'  transformed as a function of the input data into values (e.g., elements of a
#'   transition probability matrix) used for the simulation (`tparams_def` argument).
#'   The `params` argument can be omitted if the underlying parameters values are
#'  defined inside a `define_rng()` block.
#' 
#' @return `define_model()` returns an object of class `model_def`, 
#' which is a list containing the arguments to the function. `eval_model()` returns
#' a list containing [ID][id_attributes()] variables
#' identifying parameter samples, treatment strategies, patient cohorts, and time
#' intervals; the values of parameters of the transition probability matrix, 
#' utilities, and/or cost categories; the number of health states; and the number
#' of random number generation samples for the PSA. 
#' 
#' @examples 
#' 
#' # Data
#' library("data.table")
#' strategies <- data.table(strategy_id = 1:2,
#'                          strategy_name = c("Monotherapy", "Combination therapy"))
#' patients <- data.table(patient_id = 1)
#' hesim_dat <- hesim_data(strategies = strategies,
#'                        patients = patients)
#' data <- expand(hesim_dat)
#' 
#' # Model parameters
#' rng_def <- define_rng({
#'   alpha <- matrix(c(1251, 350, 116, 17,
#'                     0, 731, 512, 15,
#'                     0, 0, 1312, 437,
#'                     0, 0, 0, 469),
#'                   nrow = 4, byrow = TRUE)
#'   rownames(alpha) <- colnames(alpha) <- c("A", "B", "C", "D")
#'   lrr_mean <- log(.509)
#'   lrr_se <- (log(.710) - log(.365))/(2 * qnorm(.975))
#'   
#'   list(
#'     p_mono = dirichlet_rng(alpha),
#'     rr_comb = lognormal_rng(lrr_mean, lrr_se),
#'     u = 1,
#'     c_zido = 2278,
#'     c_lam = 2086.50,
#'     c_med = gamma_rng(mean = c(A = 2756, B = 3052, C = 9007),
#'                       sd = c(A = 2756, B = 3052, C = 9007))
#'   )
#' }, n = 2)
#'
#' tparams_def <- define_tparams({
#'   rr = ifelse(strategy_name == "Monotherapy", 1, rr_comb)
#'   list(
#'     tpmatrix = tpmatrix(
#'       C, p_mono$A_B * rr, p_mono$A_C * rr, p_mono$A_D * rr,
#'       0, C, p_mono$B_C * rr, p_mono$B_D * rr,
#'       0, 0, C, p_mono$C_D * rr,
#'       0, 0, 0, 1),
#'     utility = u,
#'     costs = list(
#'       drug = ifelse(strategy_name == "Monotherapy",
#'                     c_zido, c_zido + c_lam),
#'       medical = c_med
#'     ) 
#'   )
#' })
#' 
#' # Simulation
#' ## Define the economic model
#' model_def <- define_model(
#'   tparams_def = tparams_def,
#'   rng_def = rng_def)
#'
#' ### Evaluate the model expression to generate model inputs
#' ### This can be useful for understanding the output of a model expression
#' eval_model(model_def, data)
#' 
#' ## Create an economic model with a factory function
#' econmod <- create_CohortDtstm(model_def, data)
#'
#' @seealso [define_tparams()], [define_rng()]
#' @export
define_model <- function(tparams_def, rng_def, params = NULL,
                         n_states = NULL){
  if (!inherits(tparams_def, "list")){
    tparams_def <- list(tparams_def)
  }
  if (is.null(rng_def) & is.null(params)) {
    stop("'rng_def' and 'params' cannot both be NULL.")
  }
  x <- list(tparams_def = tparams_def, rng_def = rng_def, params = params,
            n_states = n_states)
  class(x) <- "model_def"
  check(x)
  return(x)
}

check.model_def <- function(x){
  if (!all(sapply(x$tparams_def, function (z) inherits(z, "tparams_def")))){
    stop(paste0("tparams_def must either be of class 'tparams_def'",
                "or a list of objects of class 'tparams_def'"),
         call. = FALSE)
  }
  if (!is.null(x$rng_def)) check_is_class(x$rng_def, class = "rng_def")
  if (!is.null(x$n_states)) check_scalar(x$n_states, "n_states")
}


split_params <- function(params, params_names, params_class){
  #### Split params into a list
  as_class <- function(class, z){
    switch(class,
           matrix = as.matrix(z),
           data.table = as.data.table(z),
           data.frame = as.data.frame(z),
           numeric = z[[1]],
           integer = z[[1]],
           z)
  }
  params <- split.default(params, params_names)
  for (i in 1:length(params)){
    if (ncol(params[[i]]) > 1){
      setnames(params[[i]], colnames(params[[i]]),
               gsub(pattern = paste0(names(params)[i], "."),
                    replacement = "",
                     x = colnames(params[[i]])))
    }
    params[[i]] <- as_class(params_class[[names(params)[i]]][1], params[[i]])
  }   
  return(params)
}

#' @param x An object of class `model_def` created with `define_model()`.
#' @param input_data An object of class [expanded_hesim_data][expand.hesim_data()] 
#' expanded by patients and treatment strategies. 
#' @rdname define_model
#' @export
eval_model <- function(x, input_data){
  time <- time_start <- NULL
  data <- data.table(input_data)
  
  # Step 1: RNG for parameters
  if (!is.null(x$rng_def)) {
    params <- eval_rng(x$rng_def, x$params, check = TRUE)
    n_samples <- x$rng_def$n
  } else {
    if(!attr(x$params, "checked")) check(x$params)
    params <- x$params
    n_samples <- attr(x$params, "n")
  }
  params <- as.list(params)
  params_len <- lapply(params, function (z) if (is_1d_vector(z)) 1 else ncol(z))
  params_class <- lapply(params, class)
  params_names <- rep(names(params), params_len)
  
  # Step 2: Vectorized parameter and transformed parameter values
  ## Step 2a: Parameter values
  ### Expand rows by treatment strategies and patients 
  n_obs <- nrow(data)
  
  #### First expand 'data'
  data <- data[rep(1:nrow(data), times = n_samples)]
  data <- data.table(sample = rep(1:n_samples, each = n_obs), 
                     data)
  
  #### Then expand 'params'
  params_dt <- as.data.table(params)
  params_dt <- params_dt[rep(1:nrow(params_dt), each = n_obs)]
  
  #### Get parameter values back as a list
  params <- split_params(params_dt, params_names, params_class)
  
  ## Step 2b: Transform parameters
  transform_params <- function(tparams_def){
    data_params <- cbind(data, params_dt)

    ## Expand by number of time intervals
    if (!is.null(tparams_def$times)){
      data_params <- data.table(merge(data.frame(time = tparams_def$times), 
                                      data.frame(data_params),
                                      by = NULL, sort = FALSE))
    }
    
    ## Re-split
    ### data
    data <- data_params[, c(colnames(data), "time"), with = FALSE]
    
    ### Params
    params <- data_params[, colnames(params_dt), with = FALSE]
    data_params <- NULL # Remove combined dataset
    params <- split_params(params, params_names, params_class)

    ## Evaluate parameters
    tparams <- eval_tparams(tparams_def, data, params)
    
    ## ID variables
    id_dt <- data[, c("sample", "strategy_id", "patient_id", "time"),
                  with = FALSE]
    
    ### Starting time of time intervals
    id_dt[, ("time_start") := shift(time), 
          by = c("sample", "strategy_id", "patient_id")]
    id_dt[is.na(time_start), time_start := 0]
    
    ### Other ID variables
    if (!is.null(data[["grp_id"]])) id_dt[, ("grp_id") := data[["grp_id"]]]
    if (!is.null(data[["patient_wt"]])) id_dt[, ("patient_wt") := data[["patient_wt"]]]
    
    ## Return
    return(list(id = id_dt,
                tpmatrix = tparams[["tpmatrix"]],
                utility = tparams[["utility"]],
                costs = tparams[["costs"]]))
 
  } # End function
  id <- costs <- vector(mode = "list", length = length(x$tparams_def))
  tpmatrix <- utility <- NULL
  for (i in 1:length(x$tparams_def)){
    tparams_i <- transform_params(x$tparams_def[[i]])
    id[[i]] <- tparams_i$id
    
    ### Transition probability matrix
    if (!is.null(tparams_i$tpmatrix)){
      x$n_states <- sqrt(ncol(tparams_i$tpmatrix))
      if (!is_whole_number(x$n_states)){
        stop("tpmatrix in define_tparams() must be a square matrix.", 
             call. = FALSE)
      }
      tpmatrix <- tparams_i$tpmatrix
      attr(tpmatrix, "id_index") <- i
    }
    
    ### utility
    if (!is.null(tparams_i$utility)){
      utility = tparams_i$utility
      attr(utility, "id_index") <- i
    }
    
    ### Costs
    if (!is.null(tparams_i$costs)){
      costs[[i]] <- tparams_i$costs
      for (j in 1:length(costs[[i]])){
        attr(costs[[i]][[j]], "id_index") <- i
      } 
    }
  } # End tparams_def loop
  costs <- do.call("c", costs)
  
  # Step 3: Return
  res <- list(id = id,
              tpmatrix = tpmatrix,
              utility = utility,
              costs = costs,
              n_states = x$n_states,
              n = n_samples)
  class(res) <- "eval_model" 
  check(res)
  return(res)
}

check.eval_model <- function(x){
  # Number of states
  ## Can't be NULL
  if (is.null(x$n_states)){
    stop("'n_states' cannot be NULL.", call. = FALSE)
  }
  
  ## Correct number 
  check_n_states <- function(z, name){
    if (length(dim(z)) == 2){
      if (ncol(z) != (x$n_states - 1)){
        stop(paste0("The number of columns in ", name, " must equal ",
                    "'n_states' - 1."),
             call. = FALSE)
      }
    }
  }
  check_n_states(x$utility, "'utility'")
  lapply(x$costs, check_n_states, "each element of 'costs'")
}

