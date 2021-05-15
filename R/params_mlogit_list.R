#' Parameters of a list of multinomial logit models
#' 
#' Create a list containing the parameters of multiple fitted multinomial logit models.
#' @param ... Objects of class [`params_mlogit`], which can be named.
#' 
#' @return An object of class `params_mlogit_list`, which is a list containing 
#' [`params_mlogit`] objects.
#' @export
params_mlogit_list <- function(...){
  return(check(new_params_list(..., inner_class = "params_mlogit", 
                               new_class = "params_mlogit_list")))
}

#' @export
#' @rdname create_params
create_params.multinom_list <- function(object, n = 1000, uncertainty = c("normal", "none"), ...){
  return(create_params_list(object, n = n, uncertainty = uncertainty, 
                            inner_class = "params_mlogit", new_class = "params_mlogit_list",
                            ...))
}