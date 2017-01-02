#' Posterior predictive probabilities for MCMCpack multinomial logit
#'
#' Returns sample for posterior predictive distribution from a Bayesian multinomial logit
#' fit using MCMCpack.
#'
#' @param beta Posterior distribution of parameters from Bayeisan multinomial logit model using
#' MCMCpack
#' @param newdata New data matrix of explanatory variable to make predictions with.
#' @param ncat Number of categories in multinomial logit model.
#'
#' @return List of matrices containing posterior predictive probabilities for each random sample
#' of the posterior distribution of parameters. Each matrix contains a row for each observation in
#' design matrix x and a column for each category.
#'
#' @export
#' @keywords internal
predict_MCMCmnl <- function(beta, newdata, ncat){
  return(predict_MCMCmnlC(beta, x, ncat))
}
