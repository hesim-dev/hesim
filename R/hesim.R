#' @importFrom data.table as.data.table
#' @export
data.table::as.data.table

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' @description 
#' To learn more about `hesim` visit the [website](https://hesim-dev.github.io/hesim/).
#'
#' @name hesim
#' @useDynLib hesim
#' @import data.table
#' @importFrom Rcpp evalCpp
#' @importFrom ggplot2 .data
#' @keywords internal
"_PACKAGE"