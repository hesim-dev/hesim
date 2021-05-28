#' @title Summary of a <%= class %> object
#' 
#' @description Summarize the coefficients of a parameter object by 
#' computing the point estimate, lower confidence limit, and upper confidence 
#' limit for each model term. The point estimate is the mean of the samples
#' of the coefficients and the lower and upper confidence limits are determined
#' by the `prob` argument. This is a convenient way to check whether a 
#' a parameter object has been specified correctly and sampling distributions
#' of the coefficients are as expected. 
#' 
#' @param  object An object of the appropriate class.
#' @param prob A numeric scalar in the interval `(0,1)` giving the confidence 
#' interval for coefficients. Default is 0.95 for a 95 percent interval, in which case
#' the lower and upper limits are computed using the 2.5th and 97.5th percentiles.
#' @param ... Additional arguments affecting the summary. Currently unused. 
#' 
#' @md
