
#' @title Launch an example shiny application using hesim
#'
#' @description
#' An example of how hesim can be used in a shiny application.
#'
#' @export
#'
#' @examples
#' if (interactive()) {
#'
#'  shinyIndivCtstm()
#'
#' }
shinyIndivCtstm <- function() { 
  if (!requireNamespace(package = "shiny"))
    message("Package 'shiny' is required to run this function")
  if (!requireNamespace(package = "shinydashboard"))
    message("Package 'shinydashboard' is required to run this function")
  if (!requireNamespace(package = "shinycssloaders"))
    message("Package 'shinycssloaders' is required to run this function")
  if (!requireNamespace(package = "shinyjs"))
    message("Package 'shinyjs' is required to run this function")
  if (!requireNamespace(package = "magrittr"))
    message("Package 'magrittr' is required to run this function")
  if (!requireNamespace(package = "heemod"))
    message("Package 'heemod' is required to run this function")
  if (!requireNamespace(package = "survminer"))
    message("Package 'survminer' is required to run this function")
  if (!requireNamespace(package = "DT"))
    message("Package 'DT' is required to run this function")
  if (!requireNamespace(package = "diagram"))
    message("Package 'diagram' is required to run this function")
  if (!requireNamespace(package = "rmarkdown"))
    message("Package 'rmarkdown' is required to run this function")
  if (!requireNamespace(package = "bookdown"))
    message("Package 'bookdown' is required to run this function")
  if (!requireNamespace(package = "knitr"))
    message("Package 'knitr' is required to run this function")
  if (!requireNamespace(package = "kableExtra"))
    message("Package 'kableExtra' is required to run this function")

  shiny::shinyAppDir(system.file("examples/IndivCtstm app", package = "hesim", mustWork = TRUE))
}
