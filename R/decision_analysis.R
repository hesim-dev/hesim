#' Cost-effectivenss analysis
#'
#' Conduct cost-effectiveness analysis from simulation output
#'
#' @param x Matrix containing information on mean costs and effectiveness for each simulation.
#' Should be in long form with unit of observation as the simulation and treatment arm.
#' Should have the following columns (
#' sim = simulation number,
#' arm = treatment arm,
#' c = summary of cost for each simulation and treatment arm
#' e = summary of clinicial effectiveness for each simulation and treatment arm
#' )
#' @param k Vector of willingess to pay values.
#' @param sim Name of column denoting simulation number. Default is "sim".
#' @param arm Name of column denoting treatment arm. Default is "arm".
#' @param e Name of column denoting clinical effectiveness. Default is "e".
#' @param c Name of column denoting costs. Default is "c".
#' @param custom_vars Character vector of variable names to use in creating a
#' custom summary table. Table will contain means and 95% credible intervals for
#' each variable. Can contain e and c.
#' @param custom_fun Function to apply to custom_vars to make custom table.
#' @param custom_wide. If false, arms in custom table are displayed in rows and outcomes
#'  are displayed columnwise; if TRUE, arms are dispalyed in columns and outcomes are displayed
#' rowise. Default is FALSE.
#' @return list
#'
#' @export
cea <- function(x, k, sim = "sim", arm = "arm", e = "e", c = "c",
                custom_vars = NULL, custom_fun = NULL){

  nsims <- length(unique(x[[sim]]))
  narms <- length(unique(x[[arm]]))
  setorderv(x, c(sim, arm))

  # estimates
  mce <- mce(x, k, arm, e, c, nsims, narms)
  evpi <- evpi(x, k, sim, arm, e, c, nsims, narms)
  cea.table <- cea_table(x, sim, arm, e, c, ICER = FALSE)
  l <- list(summary = cea.table, mce = mce, evpi = evpi)
  if (!is.null(custom_vars)){
    custom.table <- custom_table(x, arm, custom_vars, custom_fun)
    l <- c(l, list(custom.table = custom.table))
  }
  return(l)
}

#' Cost-effectiveness analysis with pairwise comparator
#'
#' Conduct cost-effectiveness analysis from simulation output with pairwise comparator
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
#' @param e Name of column denoting clinical effectiveness. Default is "e".
#' @param c Name of column denoting costs. Default is "c".
#' @param custom_vars Character vector of variable names to use in creating a
#' custom summary table. Table will contain means and 95% credible intervals for
#' each variable. Can contain e and c.
#' @param custom_fun Function to apply to custom_vars to make custom table.
#' @param custom_wide. If false, arms in custom table are displayed in rows and outcomes
#'  are displayed columnwise; if TRUE, arms are dispalyed in columns and outcomes are displayed
#' rowise. Default is FALSE.
#' @return list
#'
#' @export
cea_pw <- function(x, k, control, sim = "sim", arm = "arm", e = "e", c = "c",
                   custom_vars = NULL, custom_fun = NULL){

  setorderv(x, c(arm, sim))
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

  # estimates
  if (!is.null(custom_vars)){
    outcomes <- c(e, c, custom_vars[!custom_vars %in% c(e, c)])
    delta <- incr_change_calc(sim_treat, sim_control, sim, arm, outcomes, nsims, narms)
  } else{
    outcomes <- c(e, c)
    delta <- incr_change_calc(sim_treat, sim_control, sim, arm, outcomes, nsims, narms)
  }
  ceac <- ceac(sim_treat, sim_control, k, sim, arm, e, c, nsims, narms)
  einb <- einb(delta, k, sim, arm, e = paste0("i", e), c = paste0("i", c))
  cea.table <- cea_table(delta, sim, arm, e = paste0("i", e), c = paste0("i", c),
                         ICER = TRUE)
  l <- list(summary = cea.table, delta = delta, ceac = ceac, einb = einb)
  if (!is.null(custom_vars)){
    custom.table <- custom_table(delta, arm, paste0("i", custom_vars),
                                 custom_fun)
    l <- c(l, list(custom.table = custom.table))
  }
  return(l)
}

#' Incremental changes
#'
#' Calculate incremental change for outcome variables.
#'
#' @param x Matrix containing information on outcome variables for each simulation and arm.
#' @param control Name of control arm.
#' @param sim Name of column denoting simulation number.
#' @param arm Name of column denoting treatment arm.
#' @param outcomes Name of columns to calculate incremental changes for.
#' @return data.table
#'
#' @export
incr_change <- function(x, control, sim, arm, outcomes){
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
  return(incr_change_calc(x_treat, x_control, sim, arm, outcomes, nsims, narms))
}

#' Custom CEA summary table
#'
#' Custom summary table for CEA
#'
#' @param x Matrix containing information on costs and/or clinicial effectiveness for each simulation.
#' Should be in long form with unit of observation as the simulation and treatment arm.
#' @param arm Name of column denoting treatment arm. Default is "arm"
#' @param custom_vars Character vector of variable names to estimate function on.
#' @param fun List of sample statistics to estimate on each variable in custom_vars by treatment arm.
#'  Default estimates means and 95% credible intervals for each variable.
#' @return data.table data.table of summary statistics by treatment arm.
#'
#' @export
custom_table <- function(x, arm = "arm", custom_vars, FUN){
  if (is.null(FUN)){
    FUN <- function (x){
      return(list(mean = mean(x), quant = quantile(x, c(.025, .975))))
    }
  }
  ret <- dt_agg(x, FUN = FUN, by = arm, custom_vars)
  return(ret)
}

# incremental change calculation
incr_change_calc <- function(x_treat, x_control, sim, arm, outcomes, nsims, narms){
  outcomes.mat <- matrix(NA, nrow = nrow(x_treat), ncol = length(outcomes))
  colnames(outcomes.mat) <- paste0("i", outcomes)
  for (i in 1:length(outcomes)){
    outcomes.mat[, i] <- incr_changeC(x_treat[[outcomes[i]]],
                                      x_control[[outcomes[i]]],
                                      nsims, narms)
  }
  dt <- data.table(x_treat$sim, x_treat$arm, outcomes.mat)
  setnames(dt, c(sim, arm, paste0("i", outcomes)))
}

# Probability of being most cost-effective
mce <- function(x, k, arm, e, c, nsims, narms, krep, armrep){
  krep <- rep(k, each = narms)
  armrep <- rep(unique(x[[arm]]), times = length(k))
  prob.vec <- mceC(k, x[[e]], x[[c]], nsims, narms)
  prob <- data.table(krep, armrep, prob.vec)
  setnames(prob, c("k", arm, "prob"))
  return(prob)
}

# Cost effectiveness acceptability curve
ceac <- function(x_treat, x_control, k, sim, arm, e, c, nsims, narms){
  krep <- rep(k, each = narms)
  armrep <- rep(unique(x_treat[[arm]]), times = length(k))
  prob.vec <- ceacC(k, x_treat[[e]], x_treat[[c]],
                          x_control[[e]], x_control[[c]],
                          nsims, narms)
  prob <- data.table(krep, armrep, prob.vec)
  setnames(prob, c("k", arm, "prob"))
  return(prob)
}

# Expected net benefit
enb <- function(x, k, sim, arm, e, c){
  m <- x[, .(e = mean(get(e)), c = mean(get(c))), by = arm]
  enb <- matrix(NA, nrow = length(k), ncol = nrow(m))
  for (i in 1:nrow(m)){
    enb[, i] <- k * m[i, e] - m[i, c]
  }
  enb <- data.table(rep(m[[arm]], each = length(k)), rep(k, nrow(m)), c(enb))
  setnames(enb, c(arm, "k", "enb"))
  return(enb)
}

#' Expected incremental benefit
einb <- function(ix, k, sim, arm, e, c){
  einb <- enb(ix, k, sim, arm, e, c)
  setnames(einb, colnames(einb), c(arm, "k", "einb"))
  return(einb)
}

#' Expected value of perfect information
evpi <- function(x, k, sim, arm, e, c, nsims, narms){

  # Choose treatment by maximum expected benefit
  x.enb <- enb(x, k, sim, arm, e, c)
  x.enb <- dcast(x.enb, k ~ arm, value.var = "enb")[, -1, with = F]
  mu <- rowmaxC(as.matrix(x.enb))
  mu.ind <- c(rowmax_indC(as.matrix(x.enb))) + 1

  # calculate expected value of perfect information
  Vstar <- VstarC(k, x[[e]], x[[c]], nsims, narms)
  evpi <- Vstar - c(mu)
  return(data.table(k = k, evpi = evpi, ebpi = Vstar, eb = c(mu), best = mu.ind))
}

# CEA summary table
cea_table <- function(x, sim, arm, e, c, ICER = FALSE){
    FUN <- function (x){
      return(list(mean = mean(x), quant = quantile(x, c(.025, .975))))
    }
  ret <- dt_agg(x, FUN = FUN, by = arm, c(e, c))
  setnames(ret, colnames(ret), c(arm,
                                 paste0(e, c("_mean", "_lower", "_upper")),
                                 paste0(c, c("_mean", "_lower", "_upper"))))
  if (ICER == TRUE){
    ret$icer <- ret[, 5, with = F]/ret[, 2, with = F]
    ret$icer <- ifelse(ret[, 5, with = F] < 0 & ret[, 2, with = F] >= 0, "Dominates",
                       ret$icer)
    ret$icer <- ifelse(ret[, 5, with = F] < 0 & ret[, 2, with = F] < 0, "NA",
                       ret$icer)
    ret$icer <- ifelse(ret[, 5, with = F] > 0 & ret[, 2, with = F] <= 0, "Dominated",
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
