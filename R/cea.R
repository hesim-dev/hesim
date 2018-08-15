#'  A cost-effectiveness object
#'
#' An object that summarizes simulated measures of clinical effectiveness and costs from a simulation model for use in a cost-effectiveness analysis.
#'
#' 
#' @format 
#' A list containing two elements:
#' \itemize{
#' \item{costs}{ Total (discounted) costs by category.}
#' \item{QALYs}{ (Discounted) quality-adjusted life-years.}
#' }
#' 
#' @section Costs:
#' A 'costs' \code{\link{data.table}} contains the following columns:
#' \describe{
#' \item{category}{The cost category.}
#' \item{dr}{The discount rate.}
#' \item{sample}{A randomly sampled parameter set from the probabilistic sensitivity analysis (PSA)}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{grp}{An optional column denoting a subgroup. If not included, it is assumed that a single subgroup is being analyzed.}
#' \item{costs}{Costs.}
#' }
#' 
#' @section Quality-adjusted life-years:
#' A 'qalys' \code{\link{data.table}} contains the following columns:
#' \describe{
#' \item{dr}{The discount rate.}
#' \item{sample}{A randomly sampled parameter set from the probabilistic sensitivity analysis (PSA)}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{grp}{An optional column denoting a subgroup. If not included, it is assumed that a single subgroup is being analyzed.}
#' \item{qalys}{Quality-adjusted life-years}
#' }
#' 
#' @name ce
NULL

#' Individualized cost-effectiveness analysis
#'
#' Conduct Bayesian cost-effectiveness analysis (e.g. summarize a probabilistic 
#' sensitivity analysis (PSA)) by subgroup.
#' \itemize{
#'  \item \code{icea()} compute the probability that
#' each treatment is most cost-effective, the expected value of perfect
#' information, and the net monetary benefit for each treatment.
#' \item \code{icea_pw()} compares interventions to a comparator. Computed
#' quantities include the incremental cost-effectiveness ratio, the 
#' incremental net monetary benefit, output for a cost-effectiveness plane,
#' and output for a cost-effectiveness acceptability curve.
#' }
#'  
#'
#' @param x A \code{data.frame} or \code{data.table} containing simulation output with  
#' information on mean costs and effectiveness for each randomly sampled parameter set from
#' a PSA. Each row should denote a randomly sampled parameter set
#' and treatment strategy.
#' @param k Vector of willingness to pay values
#' @param comparator Name of the comparator strategy in x.
#' @param sample Name of column denoting a randomly sampled parameter set. Default is "sample".
#' @param strategy Name of column denoting treatment strategy. Default is "strategy".
#' @param grp Name of column denoting subgroup. Default is "grp".
#' @param e Name of column denoting clinical effectiveness. Default is "e".
#' @param c Name of column denoting costs. Default is "c".
#' @return \code{icea} returns a list containing four data.tables:
#' 
#' \describe{
#'   \item{summary}{A data.table of the mean, 2.5\% quantile, and 97.5\% 
#'   quantile by strategy and group for clinical effectiveness and costs.}
#'   \item{mce}{The probability that each strategy is the most effective treatment
#'   for each group for the range of specified willingness to pay values.}
#'   \item{evpi}{The expected value of perfect information by group for the range
#'   of specified willingness to pay values.}
#'    \item{nmb}{The mean, 2.5\% quantile, and 97.5\% quantile of (monetary) net benefits
#'    for the range of specified willingness to pay values.}
#' }
#' 
#' \code{icea_pw} also returns a list containing four data.tables:
#'  \describe{
#'   \item{summary}{A data.table of the mean, 2.5\% quantile, and 97.5\% 
#'   quantile by strategy and group for clinical effectiveness and costs.}
#'   \item{delta}{Incremental effectiveness and incremental cost for each simulated
#'   parameter set by strategy and group. Can be used to plot a cost-effectiveness plane. }
#'   \item{ceac}{Values needed to plot a cost-effectiveness acceptability curve by
#'   group. In other words, the probability that each strategy is more cost-effective than
#'   the comparator for the specified willingness to pay values.}
#'    \item{inmb}{The mean, 2.5\% quantile, and 97.5\% quantile of (monetary) 
#'    incremental net benefits for the range of specified willingness to pay values.}
#' }
#' @name icea
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
#' # icea
#' icea_dt <- icea(sim, k = seq(0, 200000, 500), sample = "sample", strategy = "strategy",
#'  grp = "grp", e = "e", c = "c")
#' names(icea_dt)
#' # the probability that each strategy is the most cost-effective 
#' # in each group with a willingness to pay of 20,000
#' library("data.table")
#' icea_dt$mce[k == 20000]
#' 
#' # icea_pw
#' icea_pw_dt <-  icea_pw(sim,  k = seq(0, 200000, 500), comparator = "Strategy 1",
#'                        sample = "sample", strategy = "strategy", e = "e", c = "c")
#' names(icea_pw_dt)
#' # cost-effectiveness acceptability curve
#' head(icea_pw_dt$ceac[k >= 20000])
#' @export
icea <- function(x, k, sample = "sample", strategy = "strategy", grp = "grp", e = "e", c = "c"){
  if (!is.data.table(x)){
    x <- data.table(x)
  }

  n_samples <- length(unique(x[[sample]]))
  n_strategies <- length(unique(x[[strategy]]))
  n_grps <- length(unique(x[[grp]]))
  setorderv(x, c(grp, sample, strategy))

  # estimates
  nmb <- nmb_summary(x, k, strategy, grp, e, c)
  mce <- mce(x, k, strategy, grp, e, c, n_samples, n_strategies, n_grps)
  evpi <- evpi(x, k, strategy, grp, e, c, n_samples, n_strategies, n_grps, nmb)
  summary_table <- cea_table(x, strategy, grp, e, c)
  setnames(summary_table, 
           c(paste0(e, c("_mean", "_lower", "_upper")),
            paste0(c, c("_mean", "_lower", "_upper"))),
           c(paste0("e", c("_mean", "_lower", "_upper")),
             paste0("c", c("_mean", "_lower", "_upper")))
)
  l <- list(summary = summary_table, mce = mce, evpi = evpi, nmb = nmb)
  return(l)
}

#' @export
#' @rdname icea
icea_pw <- function(x, k, comparator, sample = "sample", strategy = "strategy", 
                    grp = "grp", e = "e", c = "c"){
  if (!is.data.table(x)){
    x <- data.table(x)
  }
  setorderv(x, c(grp, strategy, sample))
  if (!comparator %in% unique(x[[strategy]])){
    stop("Chosen comparator strategy is not in x")
  }

  # treatment strategies vs comparators
  indx.comparator <- which(x[[strategy]] == comparator)
  indx_treat <- which(x[[strategy]] != comparator)
  sim_comparator <- x[indx.comparator]
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
  return(l)
}

#' Incremental treatment effect
#'
#' Calculate incremental effect of all treatment strategies on outcome variables from 
#' probabilistic sensitivity analysis relative to comparator.
#'
#' @param x A \code{data.frame} or \code{data.table} containing simulation output with  
#' information on outcome variables for each randomly sampled parameter set from
#' a PSA. Each row should denote a randomly sampled parameter set
#' and treatment strategy.
#' @param comparator Name of comparator strategy.
#' @param sample Name of column denoting a randomly sampled parameter set.
#' @param strategy Name of column denoting treatment strategy.
#' @param grp Name of column denoting subgroup.
#' @param outcomes Name of columns to calculate incremental changes for.
#' @return A data.table containing the differences in the values of each variable 
#' specified in outcomes between each treatment strategy and the 
#' comparator. It is the same output from \code{delta} generated by \code{\link{icea_pw}}.
#'
#' @examples 
#'# simulation output
#'n_samples <- 100
#'sim <- data.frame(sample = rep(seq(n_samples), 4),
#'              c = c(rlnorm(n_samples, 5, .1), rlnorm(n_samples, 5, .1),
#'                     rlnorm(n_samples, 11, .1), rlnorm(n_samples, 11, .1)),
#'              e = c(rnorm(n_samples, 8, .2), rnorm(n_samples, 8.5, .1),
#'                    rnorm(n_samples, 11, .6), rnorm(n_samples, 11.5, .6)),
#'              strategy = rep(paste0("Strategy ", seq(1, 2)),
#'                            each = n_samples * 2),
#'              grp = rep(rep(c("Group 1", "Group 2"),
#'                            each = n_samples), 2)
#')
#'# calculate incremental effect of Strategy 2 relative to Strategy 1 by group
#' ie <- incr_effect(sim, comparator = "Strategy 1", sample = "sample",
#'                         strategy = "strategy", grp = "grp", outcomes = c("c", "e"))
#' head(ie)
#' @export
incr_effect <- function(x, comparator, sample, strategy, grp, outcomes){
  x <- data.table(x)
  if (!comparator %in% unique(x[[strategy]])){
    stop("Chosen comparator strategy not in x")
  }
  indx.comparator <- which(x[[strategy]] == comparator)
  indx_treat <- which(x[[strategy]] != comparator)
  x_comparator <- x[indx.comparator]
  x_treat <- x[indx_treat]
  n_strategies <- length(unique(x_treat[[strategy]]))
  n_samples <- length(unique(x_treat[[sample]]))
  n_grps <- length(unique(x_treat[[grp]]))
  setorderv(x_treat, c(strategy, sample))
  setorderv(x_comparator, c(strategy, sample))

  # estimation
  return(calc_incr_effect(x_treat, x_comparator, sample, strategy, grp, outcomes, 
                          n_samples, n_strategies, n_grps))
}

# incremental change calculation
calc_incr_effect <- function(x_treat, x_comparator, sample, strategy, grp, outcomes,
                             n_samples, n_strategies, n_grps){
  outcomes_mat <- matrix(NA, nrow = nrow(x_treat), ncol = length(outcomes))
  colnames(outcomes_mat) <- paste0("i", outcomes)
  for (i in 1:length(outcomes)){
    outcomes_mat[, i] <- C_incr_effect(x_treat[[outcomes[i]]],
                                      x_comparator[[outcomes[i]]],
                                      n_samples, n_strategies, n_grps)
  }
  dt <- data.table(x_treat[[sample]], x_treat[[strategy]], x_treat[[grp]], outcomes_mat)
  setnames(dt, c(sample, strategy, grp, paste0("i", outcomes)))
}

# Probability of being most cost-effective
mce <- function(x, k, strategy, grp, e, c, n_samples, n_strategies, n_grps){
  k_rep <- rep(k, each = n_strategies * n_grps)
  strategy_rep <- rep(unique(x[[strategy]]), times = length(k) * n_grps)
  grp_rep <- rep(rep(unique(x[[grp]]), each = n_strategies), length(k))
  prob_vec <- C_mce(k, x[[e]], x[[c]], n_samples, n_strategies, n_grps)
  prob <- data.table(k_rep, strategy_rep, grp_rep, prob_vec)
  setnames(prob, c("k", strategy, grp, "prob"))
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
  x2 = copy(x)
  setnames(x2, c(e, c), c("e", "c"))
  m <- x2[, list(e_mean = mean(e), c_mean = mean(c),
             e_lower = stats::quantile(e, .025), 
             e_upper = stats::quantile(e, .975),
             c_lower = stats::quantile(c, .025),
             c_upper = stats::quantile(c, .975)), by = c(strategy, grp)]
  enmb <- lnmb <- unmb <- matrix(NA, nrow = length(k), ncol = nrow(m))
  for (i in 1:nrow(m)){
    enmb[, i] <- k * m[["e_mean"]][i] - m[["c_mean"]][i]
    lnmb[, i] <- k * m[["e_lower"]][i] - m[["c_lower"]][i]
    unmb[, i] <- k * m[["e_upper"]][i] - m[["c_upper"]][i]
  }
  nmb <- data.table(rep(m[[strategy]], each = length(k)), rep(m[[grp]], each = length(k)),
                    rep(k, nrow(m)), c(enmb), c(lnmb), c(unmb))
  setnames(nmb, c(strategy, grp, "k", "enmb", "lnmb", "unmb"))
  return(nmb)
}

# incremental benefit summary statistics
inmb_summary <- function(ix, k, strategy, grp, e, c){
  inmb <- nmb_summary(ix, k, strategy, grp, e, c)
  setnames(inmb, colnames(inmb), c(strategy, grp, "k", "einmb", "linmb", "uinmb"))
  return(inmb)
}

# Expected value of perfect information
evpi <- function(x, k, strategy, grp, e, c, 
                 n_samples, n_strategies, n_grps, nmb){

  # Choose treatment by maximum expected benefit
  x_nmb = copy(nmb)
  f <- stats::as.formula(paste0("k", "+", grp, "~", strategy))
  x_enmb <- dcast(x_nmb, f, value.var = "enmb")
  mu <- C_rowmax(as.matrix(x_enmb[, -c(1:2), with = FALSE]))
  mu_ind <- c(C_rowmax_index(as.matrix(x_enmb[, -c(1:2), with = FALSE]))) + 1

  # calculate expected value of perfect information
  enmbpi <- C_enmbpi(k, x[[e]], x[[c]], n_samples, n_strategies, n_grps)
  evpi <- enmbpi - c(mu)
  dt <- data.table(k = rep(k, each = n_grps),
                    grp = rep(unique(x[[grp]]), times = length(k)),
                    evpi = evpi, enmbpi = enmbpi, enmb = c(mu), best = mu_ind)
  setnames(dt, "grp", grp)
  return(dt)
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

