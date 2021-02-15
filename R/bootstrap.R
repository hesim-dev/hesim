# Bootstrapping ----------------------------------------------------------------
#' Bootstrap a statistical model
#' 
#' \code{bootstrap} is a generic function for generating bootstrap replicates of the parameters
#' of a fitted statistical model.
#' @param object A statistical model.
#' @param B Number of bootstrap replications.
#' @param max_errors Maximum number of errors that are allowed when fitting statistical models
#' during the bootstrap procedure. This argument may be useful if, for instance, the model
#' fails to converge during some bootstrap replications. Default is 0.
#' @param silent Logical indicating whether error messages should be supressed. Passed to
#' the `silent` argument of [`try()`].
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' 
#' @return Sampled values of the parameters.
#' @export
bootstrap <- function (object, B, ...) {
  UseMethod("bootstrap", object)
}

#' @name bootstrap
#' @export
bootstrap.partsurvfit <- function(object, B, max_errors = 0, silent = FALSE, ...){
  n_obs <- nrow(object$data)
  n_models <- length(object$models)
  boot <- vector(mode = "list", length = n_models)
  for (j in 1:n_models){
    boot[[j]] <- matrix(NA, nrow = B, ncol = nrow(object$models[[j]]$res.t))
    colnames(boot[[j]]) <- rownames(object$models[[j]]$res.t)
  }
  
  # Bootstrap replications
  n_errors <- 0
  i <- 1
  while(i <= B){
    index <- sample(x = 1:n_obs, size = n_obs, replace = TRUE)
    data_i <- object$data[index, ]
    fits <- vector(mode = "list", length = n_models)
    error_j <- 0
    for (j in 1:n_models){
      object$models[[j]]$call$formula <- object$models[[j]]$all.formulae[[1]]
      fit <- try(stats::update(object$models[[j]],
                    #formula = object$models[[j]]$all.formulae[[1]], 
                    # anc = object$models[[j]]$all.formulae[-1], but don't think should be needed
                    data = data_i),
                 silent = silent)
      if(inherits(fit ,"try-error")){
        n_errors <- n_errors + 1
        error_j <- error_j + 1
        if (n_errors > max_errors){
              msg <- paste0("There were ", n_errors, " errors when fitting the statistical model, ",
                  "which is greater than the maximum number of allowed errors (max_errors = ",
                  max_errors, ").")
              stop(msg)
        }
        break
      } else{
        boot[[j]][i, ] <- fit$res.t[, "est"]
      }
    }
    if (error_j == 0){
      i <- i + 1
    }
  } # end bootstrap loop
  if (n_errors > 0){
    msg <- paste0("There were ", n_errors, " errors when fitting the statistical model ",
                  "during the bootstrap procedure.")
    warning(msg)
  }
  
  # Output for params_surv_list
  params_surv_list <- vector(mode = "list", length = n_models)
  coefs <- vector(mode = "list", length = n_models)
  for (j in 1:n_models){
    n_pars <- length(object$models[[j]]$dlist$pars)
    inds_j <- flexsurvreg_inds(object$models[[j]])
    coefs[[j]] <- vector(length = n_pars, mode = "list")
    names(coefs[[j]]) <- c(object$models[[j]]$dlist$pars)
    for (k in 1:n_pars){
      coefs[[j]][[k]] <- boot[[j]][, inds_j[[k]], drop = FALSE]
      colnames(coefs[[j]][[k]])[1] <- "(Intercept)"
    }
    params_surv_list[[j]] <- new_params_surv(coefs = coefs[[j]],
                                             dist = object$models[[j]]$dlist$name,
                                             n_samples = B,
                                             aux = flexsurvreg_aux(object$models[[j]]))
  }
  names(params_surv_list) <- names(object$models)
  class(params_surv_list) <- "params_surv_list"
  return(params_surv_list)
}

# Draw parameters from multivariate normal distribution ------------------------
#' Draw parameters of statistical model from multivariate normal distribution
#' 
#' \code{normboot} is a generic function for drawing parameters from a fitted 
#' statistical model from their (asymptotic) multivariate normal distribution.
#' @param object A statistical model.
#' @param B Number of draws of the parameters.
#' @param ... Further arguments passed to or from other methods.
#' @export
#' @keywords internal
normboot <- function (object, B, ...) {
  UseMethod("normboot", object)
}

#' @export
#' @name normboot
#' @keywords internal
normboot.msm <- function(object, B = 1000, ...) {
  x <- object
  
  # As in msm:::normboot.msm simulate from a multivariate normal distribution
  ## Simulate from vector of unreplicated parameters to avoid numerical 
  ## problems with rmvnorm when lots of correlations are 1
  if (!x$foundse) stop("Asymptotic standard errors not available in fitted model.")
  sim <- MASS::mvrnorm(B, x$opt$par, x$covmat[x$paramdata$optpars, x$paramdata$optpars])
  params <- matrix(nrow = B, ncol = x$paramdata$npars)  # replicate constrained parameters.
  params[, x$paramdata$optpars] <- sim
  params[, x$paramdata$fixedpars] <- rep(x$paramdata$params[x$paramdata$fixedpars], each = B)
  params <- params[, !duplicated(abs(x$paramdata$constr)), drop = FALSE][, abs(x$paramdata$constr), drop = FALSE] *
    rep(sign(x$paramdata$constr), each = B)
  return(params)
}