#' @details Disease models may either be created from a fitted statistical 
#' model or from a parameter object. In the case of the former, `input_data`
#' is a data frame like object that is used to look for variables from
#' the statistical model that are required for simulation. In this sense,
#' `input_data` is very similar to the `newdata` argument in most [predict()]
#' methods (e.g., see [predict.lm()]). In other words, variables used in the
#' [`formula`] of the statistical model must also be in `input_data`.
#' 
#' In the case of the latter, the columns of `input_data` must be named in a 
#' manner that is consistent with the parameter object. In the typical case 
#' (e.g., with [`params_surv`] or [`params_mlogit`]), the parameter object 
#' contains coefficients from a regression model, usually stored as matrix 
#' where rows index parameter samples (i.e., for a probabilistic sensitivity
#' analysis) and columns index model terms. In such instances, there must
#' be one column from `input_data` with the same name as each model term in the 
#' coefficient matrix; that is, the columns in `input_data` are matched with
#' the columns of the coefficient matrices by name. If there are model terms
#' in the coefficient matrices that are not contained in `input_data`, then
#' an error will be thrown. 
