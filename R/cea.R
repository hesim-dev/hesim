#'  A cost-effectiveness object
#'
#' An object that summarizes simulated measures of clinical effectiveness and costs from a simulation model for use in a cost-effectiveness analysis.
#'
#' 
#' @format 
#' A list containing two elements:
#' \itemize{
#' \item{costs}{ Total (discounted) costs by category.}
#' \item{qalys}{ (Discounted) quality-adjusted life-years.}
#' }
#' 
#' @section Costs:
#' The 'costs' \code{\link{data.table}} contains the following columns:
#' \describe{
#' \item{category}{The cost category.}
#' \item{dr}{The discount rate.}
#' \item{sample}{A randomly sampled parameter set from the probabilistic sensitivity analysis (PSA)}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{grp_id}{An optional column denoting a subgroup. If not included, it is assumed that a single subgroup is being analyzed.}
#' \item{costs}{Costs.}
#' }
#' 
#' @section Quality-adjusted life-years:
#' The 'qalys' \code{\link{data.table}} contains the following columns:
#' \describe{
#' \item{dr}{The discount rate.}
#' \item{sample}{A randomly sampled parameter set from the probabilistic sensitivity analysis (PSA)}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{grp_id}{An optional column denoting a subgroup. If not included, it is assumed that a single subgroup is being analyzed.}
#' \item{qalys}{Quality-adjusted life-years}
#' }
#' 
#' @name ce
NULL

#' Cost-effectiveness analysis
#'
#' Conduct cost-effectiveness analysis (CEA) given output of an economic
#' model; that is, summarize a probabilistic sensitivity analysis (PSA), possibly
#' by subgroup.
#' \itemize{
#'  \item [cea()] computes the probability that
#' each treatment is most cost-effective, output for a cost-effectiveness acceptability frontier,
#' the expected value of perfect information, and the net monetary benefit for each treatment.
#' \item [cea_pw()] conducts pairwise CEA by comparing strategies to a comparator. Computed
#' quantities include the incremental cost-effectiveness ratio, the 
#' incremental net monetary benefit, output for a cost-effectiveness plane,
#' and output for a cost-effectiveness acceptability curve.
#' }
#'  
#'
#' @param x An object of simulation output characterizing the probability distribution
#' of clinical effectiveness and costs. If the default method is used, then `x`
#' must be a `data.frame` or `data.table` containing columns of
#' mean costs and clinical effectiveness where each row denotes a randomly sampled parameter set
#' and treatment strategy.
#' @param k Vector of willingness to pay values.
#' @param comparator Name of the comparator strategy in `x`.
#' @param sample Character name of column from `x` denoting a randomly sampled parameter set.
#' @param strategy Character name of column from `x` denoting treatment strategy.
#' @param grp Character name of column from `x` denoting subgroup. If `NULL`, then
#' it is assumed that there is only one group.
#' @param e Character name of column from `x` denoting clinical effectiveness.
#' @param c Character name of column from `x` denoting costs.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return [cea()] returns a list of four `data.table` elements.
#' 
#' \describe{
#'   \item{summary}{A `data.table` of the mean, 2.5% quantile, and 97.5% 
#'   quantile by strategy and group for clinical effectiveness and costs.}
#'   \item{mce}{The probability that each strategy is the most effective treatment
#'   for each group for the range of specified willingness to pay values. In addition,
#'   the column `best` denotes the optimal strategy (i.e., the strategy with the
#'   highest expected net monetary benefit), which can be used to plot the 
#'   cost-effectiveness acceptability frontier (CEAF).}
#'   \item{evpi}{The expected value of perfect information (EVPI) by group for the range
#'   of specified willingness to pay values. The EVPI is computed by subtracting the expected net
#'   monetary benefit given current information (i.e., the strategy with the highest
#'   expected net monetary benefit) from the expected net monetary benefit given
#'   perfect information.}
#'    \item{nmb}{The mean, 2.5% quantile, and 97.5% quantile of net monetary benefits
#'    for the range of specified willingness to pay values.}
#' }
#' 
#' \code{cea_pw} also returns a list of four `data.table` elements:
#'  \describe{
#'   \item{summary}{A data.table of the mean, 2.5% quantile, and 97.5% 
#'   quantile by strategy and group for clinical effectiveness and costs.}
#'   \item{delta}{Incremental effectiveness and incremental cost for each simulated
#'   parameter set by strategy and group. Can be used to plot a cost-effectiveness plane. }
#'   \item{ceac}{Values needed to plot a cost-effectiveness acceptability curve by
#'   group. The CEAC plots the probability that each strategy is more cost-effective than
#'   the comparator for the specified willingness to pay values.}
#'    \item{inmb}{The mean, 2.5% quantile, and 97.5% quantile of
#'    incremental net monetary benefits for the range of specified willingness to pay values.}
#' }
#' @name cea
#' @examples
#' # simulation output
#' n_samples <- 100
#' sim <- data.frame(sample = rep(seq(n_samples), 4),
#'                   c = c(rlnorm(n_samples, 5, .1), rlnorm(n_samples, 5, .1),
#'                         rlnorm(n_samples, 11, .1), rlnorm(n_samples, 11, .1)),
#'                   e = c(rnorm(n_samples, 8, .2), rnorm(n_samples, 8.5, .1),
#'                         rnorm(n_samples, 11, .6), rnorm(n_samples, 11.5, .6)),
#'                   strategy = rep(paste0("Strategy ", seq(1, 2)),
#'                                  each = n_samples * 2),
#'                   grp = rep(rep(c("Group 1", "Group 2"),
#'                             each = n_samples), 2)
#')
#'
#' # cea
#' cea <- cea(sim, k = seq(0, 200000, 500), sample = "sample", strategy = "strategy",
#'              grp = "grp", e = "e", c = "c")
#' names(cea)
#' # The probability that each strategy is the most cost-effective 
#' # in each group with a willingness to pay of 20,000
#' library("data.table")
#' cea$mce[k == 20000]
#' 
#' # cea_pw
#' cea_pw <-  cea_pw(sim,  k = seq(0, 200000, 500), comparator = "Strategy 1",
#'                     sample = "sample", strategy = "strategy", grp = "grp",
#'                      e = "e", c = "c")
#' names(cea_pw)
#' 
#' # cost-effectiveness acceptability curve
#' head(cea_pw$ceac[k >= 20000])
#' 
#' # Summarize the incremental cost-effectiveness ratio
#' format(icer(cea_pw))
#' @export
cea <- function(x, ...) {
  UseMethod("cea")
}

#' @export
#' @rdname cea
cea_pw <- function(x, ...) {
  UseMethod("cea_pw")
}

check_grp <- function(x, grp){
  if (is.null(grp)){
    grp <- "grp"
    if ("grp" %in% colnames(x)){
      x[, ("grp") := NULL]
    }
    x[, (grp) := "1"] 
  } 
  return(grp)
}

#' @export
#' @rdname cea
cea.default <- function(x, k = seq(0, 200000, 500), sample, strategy, 
                         grp = NULL, e, c, ...){
  if (!is.data.table(x)){
    x <- data.table(x)
  }
  x <- copy(x)
  grp <- check_grp(x, grp)
  n_samples <- length(unique(x[[sample]]))
  n_strategies <- length(unique(x[[strategy]]))
  n_grps <- length(unique(x[[grp]]))
  setorderv(x, c(grp, sample, strategy))

  # estimates
  nmb <- nmb_summary(x, k, strategy, grp, e, c)
  enmb_best <- enmb_best(nmb, strategy, grp)
  mce <- mce(x, k, strategy, grp, e, c, n_samples, n_strategies, n_grps, enmb_best$row)
  enmb_best[, row := NULL]
  evpi <- evpi(x, k, strategy, grp, e, c, n_samples, n_strategies, n_grps, enmb_best)
  summary_table <- cea_table(x, strategy, grp, e, c)
  setnames(summary_table, 
           c(paste0(e, c("_mean", "_lower", "_upper")),
            paste0(c, c("_mean", "_lower", "_upper"))),
           c(paste0("e", c("_mean", "_lower", "_upper")),
             paste0("c", c("_mean", "_lower", "_upper")))
)
  l <- list(summary = summary_table, mce = mce, evpi = evpi, nmb = nmb)
  class(l) <- "cea"
  attr(l, "strategy") <- strategy
  attr(l, "grp") <- grp  
  return(l)
}

#' @export
#' @rdname cea
cea_pw.default <- function(x, k = seq(0, 200000, 500), comparator, 
                            sample, strategy, 
                            grp = NULL, e, c, ...){
  if (!is.data.table(x)){
    x <- data.table(x)
  } 
  x <- copy(x)
  grp <- check_grp(x, grp)
  setorderv(x, c(grp, strategy, sample))
  if (!comparator %in% unique(x[[strategy]])){
    stop("Chosen comparator strategy is not in 'x'.",
         call. = FALSE)
  }

  # treatment strategies vs comparators
  indx_comparator <- which(x[[strategy]] == comparator)
  indx_treat <- which(x[[strategy]] != comparator)
  sim_comparator <- x[indx_comparator]
  sim_treat <- x[indx_treat]
  n_strategies <- length(unique(sim_treat[[strategy]]))
  n_samples <- length(unique(sim_treat[[sample]]))
  n_grps <- length(unique(sim_treat[[grp]]))

  # estimates
  outcomes <- c(e, c)
  delta <- calc_incr_effect(sim_treat, sim_comparator, sample, strategy, grp, outcomes, 
                              n_samples, n_strategies, n_grps)
  setnames(delta, paste0("i", e), "ie")
  setnames(delta, paste0("i", c), "ic")
  ceac <- ceac(delta, k, strategy, grp, e = "ie", c = "ic",
               n_samples, n_strategies, n_grps)
  inmb <- inmb_summary(delta, k, strategy, grp, e = "ie", c = "ic")
  summary_table <- cea_table(delta, strategy, grp, e = "ie", c = "ic", icer = TRUE)
  l <- list(summary = summary_table, delta = delta, ceac = ceac, inmb = inmb)
  class(l) <- "cea_pw"
  attr(l, "strategy") <- strategy
  attr(l, "grp") <- grp
  attr(l, "comparator") <- comparator
  if (is.factor(x[[strategy]])){
    comp_pos <- which(levels(x[[strategy]]) == comparator)
  } else {
    comp_pos <- which(sort(unique(x[[strategy]])) == comparator)
  }
  attr(l, "comparator_pos") <- comp_pos  
  return(l)
}

#' @export
#' @rdname cea
#' @param dr_qalys Discount rate for quality-adjusted life-years (QALYs).
#' @param dr_costs Discount rate for costs.
cea.ce <- function(x, k = seq(0, 200000, 500), dr_qalys, dr_costs, ...){
  category <- dr <- NULL
  sim <- cbind(x$costs[category == "total" & dr == dr_costs,
                       c("sample", "strategy_id", "grp_id", "costs")],
               x$qalys[dr == dr_qalys, "qalys", with = FALSE])
  res <- cea(sim, k = k, sample = "sample", strategy = "strategy_id",
              grp = "grp_id", e = "qalys", c = "costs")
  return(res)
}

#' @export
#' @rdname cea
cea_pw.ce <- function(x, k = seq(0, 200000, 500), comparator, dr_qalys, dr_costs, ...){
  category <- dr <- NULL
  sim <- cbind(x$costs[category == "total" & dr == dr_costs,
                       c("sample", "strategy_id", "grp_id", "costs")],
               x$qalys[dr == dr_qalys, "qalys", with = FALSE])
  res <- cea_pw(sim, k = k, comparator = comparator, sample = "sample",
                 strategy = "strategy_id", grp = "grp_id",
                 e = "qalys", c = "costs")
  return(res)
}

# Probability of being most cost-effective
mce <- function(x, k, strategy, grp, e, c, n_samples, n_strategies, n_grps,
                best_row){
  k_rep <- rep(k, each = n_strategies * n_grps)
  strategy_rep <- rep(unique(x[[strategy]]), times = length(k) * n_grps)
  grp_rep <- rep(rep(unique(x[[grp]]), each = n_strategies), length(k))
  prob_vec <- C_mce(k, x[[e]], x[[c]], n_samples, n_strategies, n_grps)
  prob <- data.table(k_rep, strategy_rep, grp_rep, prob_vec)
  setnames(prob, c("k", strategy, grp, "prob"))
  prob[, ("best") := 0]
  set(prob, best_row, "best", 1)
  setcolorder(prob, c("k", strategy, grp, "best", "prob"))
  return(prob)
}

# Cost effectiveness acceptability curve
ceac <- function(delta, k, strategy, grp, e, c, n_samples, n_strategies, n_grps){
  k_rep <- rep(k, each = n_strategies * n_grps)
  strategy_rep <- rep(unique(delta[[strategy]]), times = length(k) * n_grps)
  grp_rep <- rep(rep(unique(delta[[grp]]), each = n_strategies), length(k))
  prob_vec <- C_ceac(k, delta[[e]], delta[[c]],
                          n_samples, n_strategies, n_grps)
  prob <- data.table(k_rep, strategy_rep, grp_rep, prob_vec)
  setnames(prob, c("k", strategy, grp, "prob"))
  return(prob)
}

# net benefits summary statistics
nmb_summary <- function(x, k, strategy, grp, e, c){
  nmb <- NULL # Avoid CRAN warning for global undefined variable
  nmb_dt <- data.table(strategy = rep(x[[strategy]], times = length(k)),
                       grp = rep(x[[grp]], times = length(k)),
                       k = rep(k, each = nrow(x)),
                       e = rep(x[[e]], times = length(k)),
                       c = rep(x[[c]], times = length(k)))
  nmb_dt[, "nmb" := k * e - c]
  nmb_summary <- nmb_dt[, list("enmb" = mean(nmb),
                               "lnmb" = stats::quantile(nmb, .025),
                               "unmb" = stats::quantile(nmb, .975)),
                           by = c("strategy", "grp", "k")]
  setnames(nmb_summary, old = c("strategy", "grp"),  new = c(strategy, grp))
  return(nmb_summary)
}

# incremental benefit summary statistics
inmb_summary <- function(ix, k, strategy, grp, e, c){
  inmb <- nmb_summary(ix, k, strategy, grp, e, c)
  setnames(inmb, colnames(inmb), c(strategy, grp, "k", "einmb", "linmb", "uinmb"))
  return(inmb)
}

# Compute optimal strategy and associated ENMB
enmb_best <- function(nmb, strategy, grp){
  enmb <- NULL
  ind <- nmb[, .I[which.max(enmb)], by = c("k", grp)]$V1
  res <- nmb[ind, c(strategy, grp, "k", "enmb"), with = FALSE]
  res$row <- ind
  setnames(res, strategy, "best")
  setnames(res, "enmb", "enmb_best")
  setcolorder(res, c(grp, "k", "enmb_best", "best"))
  return(res)
}

# Expected value of perfect information
evpi <- function(x, k, strategy, grp, e, c, 
                 n_samples, n_strategies, n_grps, enmb){
  evpi <- enmbpi <- enmb_best <- NULL
  
  # calculate expected value of perfect information
  enmb$enmbpi <-  C_enmbpi(k, x[[e]], x[[c]], n_samples, n_strategies, n_grps)
  enmb[, evpi := enmbpi - enmb_best]
  setnames(enmb, "enmb_best", "enmbci")
  setcolorder(enmb, c(grp, "k", "best", "enmbci", "enmbpi", "evpi"))
  return(enmb)
}

# CEA summary table
cea_table <- function(x, strategy, grp, e, c, icer = FALSE){
  FUN <- function (x){
    return(list(mean = mean(x), quant = stats::quantile(x, c(.025, .975))))
  }
  ret <- x[, as.list(unlist(lapply(.SD, FUN))),
            by = c(strategy, grp), .SDcols = c(e, c)]
  setnames(ret, colnames(ret), c(strategy, grp,
                                 paste0(e, c("_mean", "_lower", "_upper")),
                                 paste0(c, c("_mean", "_lower", "_upper"))))
  if (icer == TRUE){
    ie_mean <- paste0(e, "_mean")
    ic_mean <- paste0(c, "_mean")
    ret$icer <- ret[, ic_mean, with = FALSE]/ret[, ie_mean, with = FALSE]
    ret$icer <- ifelse(ret[, ic_mean, with = FALSE] < 0 & ret[, ie_mean, with = FALSE] >= 0, "Dominates",
                       ret$icer)
    ret$icer <- ifelse(ret[, ic_mean, with = FALSE] > 0 & ret[, ie_mean, with = FALSE] <= 0, "Dominated",
                       ret$icer)
  }
  return(ret)
}

format_cri <- function(est, lower, upper, costs = TRUE, digits){
  if (costs){
    lower <- format_costs(lower, digits = digits)
    upper <- format_costs(upper, digits = digits)
  } else{
    lower <- format_qalys(lower, digits = digits)
    upper <- format_qalys(upper, digits = digits)
  }
  paste0(est, " (",lower, ", ", upper, ")")
}

#' Incremental cost-effectiveness ratio
#'
#' Generate a tidy table of incremental cost-effectiveness ratios (ICERs) given output from 
#' [cea_pw()] with `icer()` and format for pretty printing with `format()`.
#'
#' @param x An object of class `cea_pw` returned by [cea_pw()].
#' @param prob A numeric scalar in the interval `(0,1)` giving the confidence interval.
#' Default is 0.95 for a 95 percent interval. 
#' @param k Willingness to pay quality-adjusted life-year.
#' @return `icer()` returns a tidy `data.table` with the following columns:
#' \describe{
#' \item{strategy}{The treatment strategy.}
#' \item{grp}{The subgroup.}
#' \item{outcome}{The outcome metric.}
#' \item{mean}{The point estimate comptued as the average across the PSA samples.}
#' \item{lower}{The lower limit of the confidence interval.}
#' \item{upper}{The upper limit of the confidence interval.}
#' }
#' 
#' `format()` formats the table according to the arguments passed.
#' 
#' @seealso [`cea_pw()`]
#' @export
icer <- function(x, prob = .95, k = 50000){
  ie <- ic <- imbm <- inmb_lower <- inmb_mean <- inmb_upper <- NULL 
  
  if (!inherits(x, "cea_pw")){
    stop("'x' must be an object of class 'cea_pw'",
         call. = FALSE)
  }
  alpha <- ci_alpha(prob)
  strategy <- attributes(x)$strategy
  grp <- attributes(x)$grp 
  
  
  # Estimates for costs, QALYs, and the ICER have already been computed
  tbl <- copy(x$summary)
  
  # Compute INMB given value of k passed to function
  delta <- copy(x$delta)
  delta[, inmb := k * ie - ic]
  inmb <- delta[, list(inmb_lower = stats::quantile(inmb, alpha$lower),
                       inmb_mean = mean(inmb),
                       inmb_upper = stats::quantile(inmb, alpha$upper)),
                by = c(strategy, grp)]
  tbl <- cbind(tbl, 
               inmb[, c("inmb_lower", "inmb_mean", "inmb_upper"), with = FALSE])
  tbl <- melt(tbl, id.vars = c(strategy, grp),
              measure.vars = list(c("ie_mean", "ic_mean", "inmb_mean", "icer"),
                                  c("ie_lower", "ic_lower", "inmb_lower"),
                                  c("ie_upper", "ic_upper", "inmb_upper")),
              variable.factor = FALSE,
              variable.name = "outcome",
       value.name = c("mean", "lower", "upper"))
  tbl[, outcome := factor(outcome, levels = 1:4,
                           labels = c("Incremental QALYs",
                                      "Incremental costs",
                                      "Incremental NMB",
                                      "ICER"))]
  setorderv(tbl, c(grp, strategy))
  setattr(tbl, "class", c("icer", "data.table", "data.frame"))
  setattr(tbl, "k", k)
  setattr(tbl, "strategy", strategy)
  setattr(tbl, "grp", grp)
  return(tbl[, ])
}

#' @rdname icer
#' @param digits_qalys Number of digits to use to report QALYs.
#' @param digits_costs Number of digits to use to report costs.
#' @param pivot_from Character vector denoting a column or columns used to 
#' "widen" the data. There will be one column for each value of the variables in
#' `pivot_from`. Default is to widen so there is a column for each treatment
#' strategy. Other options are the subgroup and outcome variable (or a combination of
#' them all). Set to `NULL` if you do not want to widen the table. 
#' @export
format.icer <- function(x, digit_qalys = 2, digit_costs = 0,
                        pivot_from = "strategy") {
  y <- copy(x)
  grp <- attr(x, "grp")
  strategy <- attr(x, "strategy")
  
  # Format values
  y[, value := fcase(
    outcome %in% c("Incremental costs", "Incremental NMB"),
      format_ci(mean, lower, upper, costs = TRUE,
                digits = digit_costs),
    outcome == "Incremental QALYs", 
      format_ci(mean, lower, upper, costs = FALSE,
                digits = digit_qalys),
    outcome == "ICER", format_costs(mean, digits = digit_costs)
  )]
  y[, c("mean", "lower", "upper") := NULL]
  
  # Pivot wider
  if (!is.null(pivot_from)) {
    if (missing(pivot_from)) pivot_from <- attr(x, "strategy")
    rhs <- pivot_from
    lhs <- setdiff(c(strategy, grp, "outcome"),
                   pivot_from)
    f <- paste(paste(lhs, collapse=" + "), paste(rhs, collapse=" + "),  sep=" ~ ")
    
    y <- dcast(y, f, value.var = "value")
    
    # Sort by group and then strategy
    if (all(c(grp, strategy) %in% colnames(y))) {
      setorderv(y, c(grp, strategy))
    } else if (grp %in% colnames(y)) {
      setorderv(y, grp)
    } else if (strategy %in% colnames(y)) {
      setorderv(y, strategy)
    }
  }
  
  # Return
  return(y[, ])
}
