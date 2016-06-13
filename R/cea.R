#' cea: an R package for cost-effectiveness analysis
#'
#' It has three main goals:
#'
#' \itemize{
#' \item Identify the most important data manipulation verbs and make them
#'   easy to use from R.
#' \item Provide blazing fast performance for in-memory data by writing key
#'   pieces in C++ (using Rcpp)
#' \item Use the same interface to work with data no matter where it's stored,
#'   whether in a data frame, a data table or database.
#' }
#'
#' To learn more about dplyr, start with the vignettes:
#' \code{browseVignettes(package = "dplyr")}
#'
#' @docType package
#' @name cea
#' @useDynLib cea
#' @import data.table
#' @importFrom ggplot2 ggplot
#' @importFrom Rcpp evalCpp
NULL
