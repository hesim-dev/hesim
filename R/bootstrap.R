#' Bootstrap a statistical model
#' 
#' \code{bootstrap} is a generic function for generating bootstrap replicates of the parameters
#' of a fitted statistical model.
#' @param object A statistical model.
#' @param B Number of bootstrap replications.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' 
#' @return Sampled values of the parameters.
bootstrap <- function (object, B, ...) {
  UseMethod("bootstrap")
}

#' @export
bootstrap.partsurvfit <- function(object, B){
  n.obs <- nrow(object$data)
  n.models <- length(object$models)
  boot <- vector(mode = "list", length = n.models)
  for (j in 1:n.models){
    boot[[j]] <- matrix(NA, nrow = B, ncol = nrow(object$models[[j]]$res.t))
    colnames(boot[[j]]) <- rownames(object$models[[j]]$res.t)
  }
  
  # Bootstrap replications
  for(i in 1:B){
    index <- sample(x = 1:n.obs, size = n.obs, replace = TRUE)
    data.i <- object$data[index, ]
    fits <- vector(mode = "list", length = n.models)
    for (j in 1:n.models){
      object$models[[j]]$call$formula <- object$models[[j]]$all.formulae[[1]]
      fit <- stats::update(object$models[[j]],
                    #formula = object$models[[j]]$all.formulae[[1]], 
                    # anc = object$models[[j]]$all.formulae[-1], but don't think should be needed
                    data = data.i)
      boot[[j]][i, ] <- fit$res.t[, "est"]
    }
  } # end bootstrap loop
  
  # Output for params_surv_list
  params.surv.list <- vector(mode = "list", length = n.models)
  coefs <- vector(mode = "list", length = n.models)
  for (j in 1:n.models){
    n.pars <- length(object$models[[j]]$dlist$pars)
    inds.j <- flexsurvreg_inds(object$models[[j]])
    coefs[[j]] <- vector(length = n.pars, mode = "list")
    names(coefs[[j]]) <- c(object$models[[j]]$dlist$pars)
    for (k in 1:n.pars){
      coefs[[j]][[k]] <- boot[[j]][, inds.j[[k]], drop = FALSE]
    }
    params.surv.list[[j]] <- new_params_surv(coefs = coefs[[j]],
                                             dist = object$models[[j]]$dlist$name,
                                             n_samples = B,
                                             aux = flexsurvreg_aux(object$models[[j]]))
  }
  names(params.surv.list) <- names(object$models)
  class(params.surv.list) <- "params_surv_list"
  return(params.surv.list)
}