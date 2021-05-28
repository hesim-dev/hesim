summary_params_return <- function(extra_items) {
  all_items <- list(
    model = paste0("The name of the statistical model (e.g., a survival model for ", 
                   "particular transition or endpoint)."),
    parameter = "The name of the parameter of the survival distribution.",
    term = "The regression term.",
    estimate = paste0("The estimated value of the regression term, computed as the mean from",
                      "its probability distribution."),
    lower = "The lower limit of the confidence interval for the estimate.",
    upper = "The upper limit of the confidence interval for the estimate."
  )
  
  default_item_names <- c("term", "estimate", "lower", "upper")
  items <- all_items[names(all_items) %in% c(extra_items, default_item_names)]
  
  element_to_item <-  function(name, description) {
    paste0("\\item{", name, "}{", description, "}")
  } 
  
  describe_items_list <- lapply(seq_along(all_items), function(i) {
    element_to_item(names(all_items)[i], all_items[i])
  })
  describe_items <- paste(describe_items_list, collapse = "\n")
  
  result <- paste0("\\value{", 
                   "A \\code{\\link{data.table}} with the following columns:",
                  "\\describe{",
                  describe_items, 
                  "}}")
  result
}