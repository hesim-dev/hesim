#' Individualized cost-effectiveness analysis
#'
#' Conduct Bayesian cost-effectiveness analysis (e.g. summarize a probabilistic 
#' sensitivity analysis) by subgroup.
#'
#' @param x Matrix containing information on mean costs and effectiveness for each simulation.
#' Should be in long form with unit of observation as the simulation and treatment strategy.
#' Should have the following columns (
#' sim = simulation number,
#' strategy = treatment strategy,
#' c = summary of cost for each simulation and treatment strategy
#' e = summary of clinicial effectiveness for each simulation and treatment strategy)
#' @param k Vector of willingess to pay values
#' @param comparator Name of the comparator strategy in x.
#' @param sim Name of column denoting simulation number. Default is "sim".
#' @param strategy Name of column denoting treatment strategy. Default is "strategy".
#' @param grp Name of column denoting subgroup. Default is "grp".
#' @param e Name of column denoting clinical effectiveness. Default is "e".
#' @param c Name of column denoting costs. Default is "c".
#' @param custom_vars Character vector of variable names to use in creating a
#' custom summary table. Table will contain means and 95\% credible intervals for
#' each variable. Can contain e and c.
#' @param custom_fun Function to apply to custom_vars to make custom table. If 
#' \code{custom_vars} is not NULL and \code{custom_fun} is NULL, then returns the mean,
#' 2.5\% quantile, and 97.5\% quantile for each variable in \code{custom_vars}.
#' @return \code{icea} returns a list containing four data.tables:
#' 
#' \describe{
#'   \item{summary}{A data.table of the mean, 2.5\% quantile, and 97.5\% 
#'   quantile by strategy and group for clinical effectiveness and costs.}
#'   \item{mce}{The probability that each strategy is the most effective treatment
#'   for each group for the range of specified willingess to pay values.}
#'   \item{evpi}{The expected value of perfect information by group for the range
#'   of specified willingess to pay values.}
#'    \item{nmb}{The mean, 2.5\% quantile, and 97.5\% quantile of (monetary) net benefits
#'    for the range of specified willingess to pay values.}
#' }
#' In addition, if \code{custom_vars} is not NULL, \code{icea} returns \code{custom.table}, which is
#'  a data.table containing summary statistics for each variable in \code{custom_vars}
#'   by strategy and group.
#' 
#' \code{icea_pw} also returns a list containing four data.tables:
#'  \describe{
#'   \item{summary}{A data.table of the mean, 2.5\% quantile, and 97.5\% 
#'   quantile by strategy and group for clinical effectiveness and costs.}
#'   \item{delta}{Incremental effectiveness and incremental cost for each simulated
#'   parameter set by strategy and group. Can be used to plot a cost-effectiveness plane. 
#'   Also returns the difference between each treatment strategy and the comparator for each 
#'   variable in \code{custom_vars} if \code{custom_vars} is not NULL.}
#'   \item{ceac}{Values needed to plot a cost-effectiveness acceptability curve by
#'   group. In other words, the probability that each strategy is more cost-effective than
#'   the comparator for the specified willingess to pay values.}
#'    \item{inmb}{The mean, 2.5\% quantile, and 97.5\% quantile of (monetary) 
#'    incremental net benefits for the range of specified willingess to pay values.}
#' }
#' If \code{custom_vars} is not NULL, also returns \code{custom.table}, which is
#'  a data.table containing summary statistics for the values of each variable returned
#'   in \code{delta} by strategy and group.
#' @name icea
#' @export
icea <- function(x, k, sim = "sim", strategy = "strategy", grp = "grp", e = "e", c = "c",
                custom_vars = NULL, custom_fun = NULL){
  if (!is.data.table(x)){
    x <- data.table(x)
  }

  nsims <- length(unique(x[[sim]]))
  nstrategies <- length(unique(x[[strategy]]))
  ngrps <- length(unique(x[[grp]]))
  setorderv(x, c(grp, sim, strategy))

  # estimates
  nmb <- nmb_summary(x, k, sim, strategy, grp, e, c)
  mce <- mce(x, k, strategy, grp, e, c, nsims, nstrategies, ngrps)
  evpi <- evpi(x, k, sim, strategy, grp, e, c, nsims, nstrategies, ngrps, nmb)
  cea.table <- cea_table(x, sim, strategy, grp, e, c, ICER = FALSE)
  setnames(cea.table, 
           c(paste0(e, c("_mean", "_lower", "_upper")),
            paste0(c, c("_mean", "_lower", "_upper"))),
           c(paste0("e", c("_mean", "_lower", "_upper")),
             paste0("c", c("_mean", "_lower", "_upper")))
)
  l <- list(summary = cea.table, mce = mce, evpi = evpi, nmb = nmb)
  if (!is.null(custom_vars)){
    custom.table <- custom_table(x, strategy, grp, custom_vars, custom_fun)
    l <- c(l, list(custom.table = custom.table))
  }
  return(l)
}

#' @export
#' @rdname icea
icea_pw <- function(x, k, comparator, sim = "sim", strategy = "strategy", grp = "grp", e = "e", c = "c",
                   custom_vars = NULL, custom_fun = NULL){
  if (!is.data.table(x)){
    x <- data.table(x)
  }
  setorderv(x, c(grp, strategy, sim))
  if (!comparator %in% unique(x[[strategy]])){
    stop("Chosen comparator strategy is not in x")
  }

  # treatment strategies vs comparators
  indx.comparator <- which(x[[strategy]] == comparator)
  indx.treat <- which(x[[strategy]] != comparator)
  sim_comparator <- x[indx.comparator]
  sim_treat <- x[indx.treat]
  nstrategies <- length(unique(sim_treat[[strategy]]))
  nsims <- length(unique(sim_treat[[sim]]))
  ngrps <- length(unique(sim_treat[[grp]]))

  # estimates
  if (!is.null(custom_vars)){
    outcomes <- c(e, c, custom_vars[!custom_vars %in% c(e, c)])
    delta <- calc_incr_effect(sim_treat, sim_comparator, sim, strategy, grp, outcomes, 
                              nsims, nstrategies, ngrps)
  } else{
    outcomes <- c(e, c)
    delta <- calc_incr_effect(sim_treat, sim_comparator, sim, strategy, grp, outcomes, 
                              nsims, nstrategies, ngrps)
  }
  setnames(delta, paste0("i", e), "ie")
  setnames(delta, paste0("i", c), "ic")
  ceac <- ceac(delta, k, sim, strategy, grp, e = "ie", c = "ic",
               nsims, nstrategies, ngrps)
  inmb <- inmb_summary(delta, k, sim, strategy, grp, e = "ie", c = "ic")
  cea.table <- cea_table(delta, sim, strategy, grp, e = "ie", c = "ic",
                         ICER = TRUE)
  l <- list(summary = cea.table, delta = delta, ceac = ceac, inmb = inmb)
  if (!is.null(custom_vars)){
    if (any(custom_vars == e)){
        custom_vars[custom_vars == e] <- "e"
    }
    if (any(custom_vars == c)){
      custom_vars[custom_vars == c] <- "c"
    }
    custom.table <- custom_table(delta, strategy, grp, paste0("i", custom_vars),
                                 custom_fun)
    l <- c(l, list(custom.table = custom.table))
  }
  return(l)
}

#' Incremental treatment effect
#'
#' Calculate incremental effect of all treatment strategies on outcome variables from 
#' probabilistic sensitivity analysis relative to comparator.
#'
#' @param x Matrix containing information on outcome variables for each simulation and strategy.
#' @param comparator Name of comparator strategy.
#' @param sim Name of column denoting simulation number.
#' @param strategy Name of column denoting treatment strategy.
#' @param grp Name of column denoting subgroup.
#' @param outcomes Name of columns to calculate incremental changes for.
#' @return A data.table containing the differences in the values of each variable 
#' specified in outcomes between each treatment strategy and the 
#' comparator. It is the same output generated
#' in \code{delta} from \code{icea_pw}.
#'
#' @export
incr_effect <- function(x, comparator, sim, strategy, grp, outcomes){
  if (!comparator %in% unique(x[[strategy]])){
    stop("Chosen comparator strategy not in x")
  }
  indx.comparator <- which(x[[strategy]] == comparator)
  indx.treat <- which(x[[strategy]] != comparator)
  x_comparator <- x[indx.comparator]
  x_treat <- x[indx.treat]
  nstrategies <- length(unique(x_treat[[strategy]]))
  nsims <- length(unique(x_treat[[sim]]))
  ngrps <- length(unique(x_treat[[grp]]))
  setorderv(x_treat, c(strategy, sim))
  setorderv(x_comparator, c(strategy, sim))

  # estimation
  return(calc_incr_effect(x_treat, x_comparator, sim, strategy, grp, outcomes, 
                          nsims, nstrategies, ngrps))
}


#' Custom CEA summary table
#'
#' Custom table summarizing outcomes from probabilistic sensitivity analysis.
#'
#' @param x Matrix containing information on outcome variables for each simulation and strategy.
#' @param strategy Name of column denoting treatment strategy.
#' @param grp Name of columne denoting subgroup
#' @param custom_vars Name of custom variables to calculate summary statistic for.
#' @param FUN summary statistic function.
#' @return A data.table of summary statistics for each variable specified in 
#' \code{custom_vars}. By default, returns the mean, 2.5\%, and 97.5\% quantile of
#' each variable. Different summary statistics can be calculated using FUN. 
#' This function is used in \code{icea} and \code{icea_pw} to create the
#'  \code{custom.table} output.
#'
#' @export
custom_table <- function(x, strategy, grp, custom_vars, FUN = NULL){
  if (is.null(FUN)){
    FUN <- function (x){
      return(list(mean = mean(x), quant = quantile(x, c(.025, .975))))
    }
  }
  ret <- dt_agg(x, FUN = FUN, by = c(strategy, grp), custom_vars)
  return(ret)
}

# incremental change calculation
calc_incr_effect <- function(x_treat, x_comparator, sim, strategy, grp, outcomes,
                             nsims, nstrategies, ngrps){
  outcomes.mat <- matrix(NA, nrow = nrow(x_treat), ncol = length(outcomes))
  colnames(outcomes.mat) <- paste0("i", outcomes)
  for (i in 1:length(outcomes)){
    outcomes.mat[, i] <- incr_effectC(x_treat[[outcomes[i]]],
                                      x_comparator[[outcomes[i]]],
                                      nsims, nstrategies, ngrps)
  }
  dt <- data.table(x_treat[[sim]], x_treat[[strategy]], x_treat[[grp]], outcomes.mat)
  setnames(dt, c(sim, strategy, grp, paste0("i", outcomes)))
}

# Probability of being most cost-effective
mce <- function(x, k, strategy, grp, e, c, nsims, nstrategies, ngrps){
  krep <- rep(k, each = nstrategies * ngrps)
  strategyrep <- rep(unique(x[[strategy]]), times = length(k) * ngrps)
  grprep <- rep(rep(unique(x[[grp]]), each = nstrategies), length(k))
  prob.vec <- mceC(k, x[[e]], x[[c]], nsims, nstrategies, ngrps)
  prob <- data.table(krep, strategyrep, grprep, prob.vec)
  setnames(prob, c("k", strategy, grp, "prob"))
  return(prob)
}

# Cost effectiveness acceptability curve
ceac <- function(delta, k, sim, strategy, grp, e, c, nsims, nstrategies, ngrps){
  krep <- rep(k, each = nstrategies * ngrps)
  strategyrep <- rep(unique(delta[[strategy]]), times = length(k) * ngrps)
  grprep <- rep(rep(unique(delta[[grp]]), each = nstrategies), length(k))
  prob.vec <- ceacC(k, delta[[e]], delta[[c]],
                          nsims, nstrategies, ngrps)
  prob <- data.table(krep, strategyrep, grprep, prob.vec)
  setnames(prob, c("k", strategy, grp, "prob"))
  return(prob)
}

# net benefits summary statistics
nmb_summary <- function(x, k, sim, strategy, grp, e, c){
  m <- x[, .(e_mean = mean(get(e)), c_mean = mean(get(c)),
             e_lower = quantile(get(e), .025), 
             e_upper = quantile(get(e), .975),
             c_lower = quantile(get(c), .025),
             c_upper = quantile(get(c), .975)), by = c(strategy, grp)]
  enmb <- lnmb <- unmb <- matrix(NA, nrow = length(k), ncol = nrow(m))
  for (i in 1:nrow(m)){
    enmb[, i] <- k * m[i, e_mean] - m[i, c_mean]
    lnmb[, i] <- k * m[i, e_lower] - m[i, c_lower]
    unmb[, i] <- k * m[i, e_upper] - m[i, c_upper]
  }
  nmb <- data.table(rep(m[[strategy]], each = length(k)), rep(m[[grp]], each = length(k)),
                    rep(k, nrow(m)), c(enmb), c(lnmb), c(unmb))
  setnames(nmb, c(strategy, grp, "k", "enmb", "lnmb", "unmb"))
  return(nmb)
}

# incremental benefit summary statistics
inmb_summary <- function(ix, k, sim, strategy, grp, e, c){
  inmb <- nmb_summary(ix, k, sim, strategy, grp, e, c)
  setnames(inmb, colnames(inmb), c(strategy, grp, "k", "einmb", "linmb", "uinmb"))
  return(inmb)
}

# Expected value of perfect information
evpi <- function(x, k, sim, strategy, grp, e, c, nsims, nstrategies, ngrps, nmb){

  # Choose treatment by maximum expected benefit
  x.nmb = copy(nmb)
  f <- as.formula(paste0("k", "+", grp, "~", strategy))
  x.enmb <- dcast(x.nmb, f, value.var = "enmb")
  mu <- rowmaxC(as.matrix(x.enmb[, -c(1:2), with = FALSE]))
  mu.ind <- c(rowmax_indC(as.matrix(x.enmb[, -c(1:2), with = FALSE]))) + 1

  # calculate expected value of perfect information
  Vstar <- VstarC(k, x[[e]], x[[c]], nsims, nstrategies, ngrps)
  evpi <- Vstar - c(mu)
  dt <- data.table(k = rep(k, each = ngrps),
                    grp = rep(unique(x[[grp]]), times = length(k)),
                    evpi = evpi, enmbpi = Vstar, enmb = c(mu), best = mu.ind)
  setnames(dt, "grp", grp)
  return(dt)
}

# CEA summary table
cea_table <- function(x, sim, strategy, grp, e, c, ICER = FALSE){
    FUN <- function (x){
      return(list(mean = mean(x), quant = quantile(x, c(.025, .975))))
    }
  ret <- dt_agg(x, FUN = FUN, by = c(strategy, grp), c(e, c))
  setnames(ret, colnames(ret), c(strategy, grp,
                                 paste0(e, c("_mean", "_lower", "_upper")),
                                 paste0(c, c("_mean", "_lower", "_upper"))))
  ie_mean <- paste0(e, "_mean")
  ic_mean <- paste0(c, "_mean")
  if (ICER == TRUE){
    ret$icer <- ret[, ic_mean, with = FALSE]/ret[, ie_mean, with = FALSE]
    ret$icer <- ifelse(ret[, ic_mean, with = FALSE] < 0 & ret[, ie_mean, with = FALSE] >= 0, "Dominates",
                       ret$icer)
    ret$icer <- ifelse(ret[, ic_mean, with = FALSE] > 0 & ret[, ie_mean, with = FALSE] <= 0, "Dominated",
                       ret$icer)
  }
  return(ret)
}

# Helper functions
dt_agg <- function(x, FUN, by, sdcols){
  y <- x[, as.list(unlist(lapply(.SD, FUN))),
            by = by, .SDcols = sdcols]
  return(y)
}

EVAL = function(...) {
  eval(parse(text=paste0(...)),envir=parent.frame(2))
}
