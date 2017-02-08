#' Personalized cost-effectiveness analysis
#'
#' Conduct Bayesian cost-effectiveness analysis (e.g. summarize a probabilistic 
#' sensitivity analysis) by subgroup.
#'
#' @param x Matrix containing information on mean costs and effectiveness for each simulation.
#' Should be in long form with unit of observation as the simulation and treatment arm.
#' Should have the following columns (
#' sim = simulation number,
#' arm = treatment arm,
#' c = summary of cost for each simulation and treatment arm
#' e = summary of clinicial effectiveness for each simulation and treatment arm)
#' @param k Vector of willingess to pay values
#' @param control Name of the control arm in x.
#' @param sim Name of column denoting simulation number. Default is "sim".
#' @param arm Name of column denoting treatment arm. Default is "arm".
#' @param grp Name of column denoting subgroup. Default is "grp".
#' @param e Name of column denoting clinical effectiveness. Default is "e".
#' @param c Name of column denoting costs. Default is "c".
#' @param custom_vars Character vector of variable names to use in creating a
#' custom summary table. Table will contain means and 95\% credible intervals for
#' each variable. Can contain e and c.
#' @param custom_fun Function to apply to custom_vars to make custom table. If 
#' \code{custom_vars} is not NULL and \code{custom_fun} is NULL, then returns the mean,
#' 2.5\% quantile, and 97.5\% quantile for each variable in \code{custom_vars}.
#' @return \code{pcea} returns a list containing four data.tables:
#' 
#' \describe{
#'   \item{summary}{A data.table of the mean, 2.5\% quantile, and 97.5\% 
#'   quantile by arm and group for clinical effectiveness and costs.}
#'   \item{mce}{The probability that each arm is the most effective treatment
#'   for each group for the range of specified willingess to pay values.}
#'   \item{evpi}{The expected value of perfect information by group for the range
#'   of specified willingess to pay values.}
#'    \item{nmb}{The mean, 2.5\% quantile, and 97.5\% quantile of (monetary) net benefits
#'    for the range of specified willingess to pay values.}
#' }
#' In addition, if \code{custom_vars} is not NULL, \code{pcea} returns \code{custom.table}, which is
#'  a data.table containing summary statistics for each variable in \code{custom_vars}
#'   by arm and group.
#' 
#' \code{pcea_pw} also returns a list containing four data.tables:
#'  \describe{
#'   \item{summary}{A data.table of the mean, 2.5\% quantile, and 97.5\% 
#'   quantile by arm and group for clinical effectiveness and costs.}
#'   \item{delta}{Incremental effectiveness and incremental cost for each simulated
#'   parameter set by arm and group. Can be used to plot a cost-effectiveness plane. 
#'   Also returns the difference between each treatment arm and the comparator for each 
#'   variable in \code{custom_vars} if \code{custom_vars} is not NULL.}
#'   \item{ceac}{Values needed to plot a cost-effectiveness acceptability curve by
#'   group. In other words, the probability that each arm is more cost-effective than
#'   the comparator for the specified willingess to pay values.}
#'    \item{inmb}{The mean, 2.5\% quantile, and 97.5\% quantile of (monetary) 
#'    incremental net benefits for the range of specified willingess to pay values.}
#' }
#' If \code{custom_vars} is not NULL, also returns \code{custom.table}, which is
#'  a data.table containing summary statistics for the values of each variable returned
#'   in \code{delta} by arm and group.
#' @name pcea
#' @export
pcea <- function(x, k, sim = "sim", arm = "arm", grp = "grp", e = "e", c = "c",
                custom_vars = NULL, custom_fun = NULL){
  if (!is.data.table(x)){
    x <- data.table(x)
  }

  nsims <- length(unique(x[[sim]]))
  narms <- length(unique(x[[arm]]))
  ngrps <- length(unique(x[[grp]]))
  setorderv(x, c(grp, sim, arm))

  # estimates
  nmb <- nmb_summary(x, k, sim, arm, grp, e, c)
  mce <- mce(x, k, arm, grp, e, c, nsims, narms, ngrps)
  evpi <- evpi(x, k, sim, arm, grp, e, c, nsims, narms, ngrps, nmb)
  cea.table <- cea_table(x, sim, arm, grp, e, c, ICER = FALSE)
  l <- list(summary = cea.table, mce = mce, evpi = evpi, nmb = nmb)
  if (!is.null(custom_vars)){
    custom.table <- custom_table(x, arm, grp, custom_vars, custom_fun)
    l <- c(l, list(custom.table = custom.table))
  }
  return(l)
}

#' @export
#' @rdname pcea
pcea_pw <- function(x, k, control, sim = "sim", arm = "arm", grp = "grp", e = "e", c = "c",
                   custom_vars = NULL, custom_fun = NULL){
  if (!is.data.table(x)){
    x <- data.table(x)
  }
  setorderv(x, c(grp, arm, sim))
  if (!control %in% unique(x[[arm]])){
    stop("Chosen control arm is not in x")
  }

  # treatment arms vs comparators
  indx.control <- which(x[[arm]] == control)
  indx.treat <- which(x[[arm]] != control)
  sim_control <- x[indx.control]
  sim_treat <- x[indx.treat]
  narms <- length(unique(sim_treat[[arm]]))
  nsims <- length(unique(sim_treat[[sim]]))
  ngrps <- length(unique(sim_treat[[grp]]))

  # estimates
  if (!is.null(custom_vars)){
    outcomes <- c(e, c, custom_vars[!custom_vars %in% c(e, c)])
    delta <- calc_incr_effect(sim_treat, sim_control, sim, arm, grp, outcomes, nsims, narms, ngrps)
  } else{
    outcomes <- c(e, c)
    delta <- calc_incr_effect(sim_treat, sim_control, sim, arm, grp, outcomes, nsims, narms, ngrps)
  }
  ceac <- ceac(delta, k, sim, arm, grp, e = paste0("i", e), c = paste0("i", c),
               nsims, narms, ngrps)
  inmb <- inmb_summary(delta, k, sim, arm, grp, e = paste0("i", e), c = paste0("i", c))
  cea.table <- cea_table(delta, sim, arm, grp, e = paste0("i", e), c = paste0("i", c),
                         ICER = TRUE)
  l <- list(summary = cea.table, delta = delta, ceac = ceac, inmb = inmb)
  if (!is.null(custom_vars)){
    custom.table <- custom_table(delta, arm, grp, paste0("i", custom_vars),
                                 custom_fun)
    l <- c(l, list(custom.table = custom.table))
  }
  return(l)
}

#' Incremental treatment effect
#'
#' Calculate incremental effect of all treatment arms on outcome variables from 
#' probabilistic sensitivity analysis relative to comparator.
#'
#' @param x Matrix containing information on outcome variables for each simulation and arm.
#' @param control Name of control arm (i.e. comparator).
#' @param sim Name of column denoting simulation number.
#' @param arm Name of column denoting treatment arm.
#' @param outcomes Name of columns to calculate incremental changes for.
#' @return A data.table containing the differences in the values of each variable 
#' specified in outcomes between each treatment arm and the 
#' comparator. It is the same output generated
#' in \code{delta} from \code{pcea_pw}.
#'
#' @export
incr_effect <- function(x, control, sim, arm, outcomes){
  if (!control %in% unique(x[[arm]])){
    stop("Chosen control arm is not in x")
  }
  indx.control <- which(x[[arm]] == control)
  indx.treat <- which(x[[arm]] != control)
  x_control <- x[indx.control]
  x_treat <- x[indx.treat]
  narms <- length(unique(x_treat[[arm]]))
  nsims <- length(unique(x_treat[[sim]]))
  setorderv(x_treat, c(arm, sim))
  setorderv(x_control, c(arm, sim))

  # estimation
  return(calc_incr_effect(x_treat, x_control, sim, arm, outcomes, nsims, narms))
}


#' Custom CEA summary table
#'
#' Custom table summarizing outcomes from probabilistic sensitivity analysis.
#'
#' @param x Matrix containing information on outcome variables for each simulation and arm.
#' @param arm Name of column denoting treatment arm.
#' @param grp Name of columne denoting subgroup
#' @param custom_vars Name of custom variables to calculate summary statistic for.
#' @param FUN summary statistic function.
#' @return A data.table of summary statistics for each variable specified in 
#' \code{custom_vars}. By default, returns the mean, 2.5\%, and 97.5\% quantile of
#' each variable. Different summary statistics can be calculated using FUN. 
#' This function is used in \code{pcea} and \code{pcea_pw} to create the
#'  \code{custom.table} output.
#'
#' @export
custom_table <- function(x, arm, grp, custom_vars, FUN = NULL){
  if (is.null(FUN)){
    FUN <- function (x){
      return(list(mean = mean(x), quant = quantile(x, c(.025, .975))))
    }
  }
  ret <- dt_agg(x, FUN = FUN, by = c(arm, grp), custom_vars)
  return(ret)
}

# incremental change calculation
calc_incr_effect <- function(x_treat, x_control, sim, arm, grp, outcomes,
                             nsims, narms, ngrps){
  outcomes.mat <- matrix(NA, nrow = nrow(x_treat), ncol = length(outcomes))
  colnames(outcomes.mat) <- paste0("i", outcomes)
  for (i in 1:length(outcomes)){
    outcomes.mat[, i] <- incr_effectC(x_treat[[outcomes[i]]],
                                      x_control[[outcomes[i]]],
                                      nsims, narms, ngrps)
  }
  dt <- data.table(x_treat$sim, x_treat$arm, x_treat$grp, outcomes.mat)
  setnames(dt, c(sim, arm, grp, paste0("i", outcomes)))
}

# Probability of being most cost-effective
mce <- function(x, k, arm, grp, e, c, nsims, narms, ngrps){
  krep <- rep(k, each = narms * ngrps)
  armrep <- rep(unique(x[[arm]]), times = length(k) * ngrps)
  grprep <- rep(rep(unique(x[[grp]]), each = narms), length(k))
  prob.vec <- mceC(k, x[[e]], x[[c]], nsims, narms, ngrps)
  prob <- data.table(krep, armrep, grprep, prob.vec)
  setnames(prob, c("k", arm, grp, "prob"))
  return(prob)
}

# Cost effectiveness acceptability curve
ceac <- function(delta, k, sim, arm, grp, e, c, nsims, narms, ngrps){
  krep <- rep(k, each = narms * ngrps)
  armrep <- rep(unique(delta[[arm]]), times = length(k) * ngrps)
  grprep <- rep(rep(unique(delta[[grp]]), each = narms), length(k))
  prob.vec <- ceacC(k, delta[[e]], delta[[c]],
                          nsims, narms, ngrps)
  prob <- data.table(krep, armrep, grprep, prob.vec)
  setnames(prob, c("k", arm, grp, "prob"))
  return(prob)
}

# net benefits summary statistics
nmb_summary <- function(x, k, sim, arm, grp, e, c){
  m <- x[, .(e_mean = mean(get(e)), c_mean = mean(get(c)),
             e_lower = quantile(get(e), .025), 
             e_upper = quantile(get(e), .975),
             c_lower = quantile(get(c), .025),
             c_upper = quantile(get(c), .975)), by = c(arm, grp)]
  enmb <- lnmb <- unmb <- matrix(NA, nrow = length(k), ncol = nrow(m))
  for (i in 1:nrow(m)){
    enmb[, i] <- k * m[i, e_mean] - m[i, c_mean]
    lnmb[, i] <- k * m[i, e_lower] - m[i, c_lower]
    unmb[, i] <- k * m[i, e_upper] - m[i, c_upper]
  }
  nmb <- data.table(rep(m[[arm]], each = length(k)), rep(m[[grp]], each = length(k)),
                    rep(k, nrow(m)), c(enmb), c(lnmb), c(unmb))
  setnames(nmb, c(arm, grp, "k", "enmb", "lnmb", "unmb"))
  return(nmb)
}

# incremental benefit summary statistics
inmb_summary <- function(ix, k, sim, arm, grp, e, c){
  inmb <- nmb_summary(ix, k, sim, arm, grp, e, c)
  setnames(inmb, colnames(inmb), c(arm, grp, "k", "einmb", "linmb", "uinmb"))
  return(inmb)
}

# Expected value of perfect information
evpi <- function(x, k, sim, arm, grp, e, c, nsims, narms, ngrps, nmb){

  # Choose treatment by maximum expected benefit
  x.nmb = copy(nmb)
  x.enmb <- dcast(x.nmb, k + grp ~ arm, value.var = "enmb")
  mu <- rowmaxC(as.matrix(x.enmb[, -c(1:2), with = FALSE]))
  mu.ind <- c(rowmax_indC(as.matrix(x.enmb[, -c(1:2), with = FALSE]))) + 1

  # calculate expected value of perfect information
  Vstar <- VstarC(k, x[[e]], x[[c]], nsims, narms, ngrps)
  evpi <- Vstar - c(mu)
  return(data.table(k = rep(k, each = ngrps),
                    grp = rep(unique(x[[grp]]), times = length(k)),
                    evpi = evpi, enmbpi = Vstar, enmb = c(mu), best = mu.ind))
}

# CEA summary table
cea_table <- function(x, sim, arm, grp, e, c, ICER = FALSE){
    FUN <- function (x){
      return(list(mean = mean(x), quant = quantile(x, c(.025, .975))))
    }
  ret <- dt_agg(x, FUN = FUN, by = c(arm, grp), c(e, c))
  setnames(ret, colnames(ret), c(arm, grp,
                                 paste0(e, c("_mean", "_lower", "_upper")),
                                 paste0(c, c("_mean", "_lower", "_upper"))))
  if (ICER == TRUE){
    ret$icer <- ret[, icost_mean]/ret[, iqalys_mean]
    ret$icer <- ifelse(ret[, icost_mean] < 0 & ret[, iqalys_mean] >= 0, "Dominates",
                       ret$icer)
    ret$icer <- ifelse(ret[, icost_mean] < 0 & ret[, iqalys_mean] < 0, "NA",
                       ret$icer)
    ret$icer <- ifelse(ret[, icost_mean] > 0 & ret[, iqalys_mean] <= 0, "Dominated",
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
