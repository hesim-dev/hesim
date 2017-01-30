#' US Male Lifetables 2011
#'
#' Life tables by single-year of age from National Vital Statistics Reports
#' Volume 64, Number 11.
#'
#' @name lifetable
#' @format A data frame with 101 rows and 3 variables:
#' \describe{
#'   \item{age}{age in years}
#'   \item{qx}{Probability of dying between ages x and x + 1}
#'   \item{lx}{Number surviving to age x}
#'   \item{dx}{Number dying between ages x and x + 1}
#'   \item{L}{Person-years lived between ages x and x + 1}
#'   \item{Tx}{Total number of person-years lived above age x}
#'   \item{ex}{Expectation of life at age x}
#'
#' }
#' @source \url{http://www.cdc.gov/nchs/data/nvsr/nvsr64/nvsr64_11.pdf}
NULL

#' @rdname lifetable
"lifetable_male"

#' @rdname lifetable
"lifetable_female"
