# Incremental treatment effect -------------------------------------------------
#' Incremental treatment effect
#'
#' Computes incremental effect for all treatment strategies 
#' on outcome variables from a probabilistic sensitivity analysis relative to a comparator.
#'
#' @param x A \code{data.frame} or \code{data.table} containing simulation output with  
#' information on outcome variables for each randomly sampled parameter set from
#' a PSA. Each row should denote a randomly sampled parameter set
#' and treatment strategy.
#' @param comparator The comparator strategy. If the strategy column is a character
#' variable, then must be a character; if the strategy column is an integer variable,
#' then must be an integer.
#' @param sample Character name of column denoting a randomly sampled parameter set.
#' @param strategy Character name of column denoting treatment strategy.
#' @param grp Character name of column denoting subgroup. If \code{NULL}, then
#' it is assumed that there is only one group.
#' @param outcomes Name of columns to compute incremental changes for.
#' @return A \code{data.table} containing the differences in the values of each variable 
#' specified in outcomes between each treatment strategy and the 
#' comparator. 
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
incr_effect <- function(x, comparator, sample, strategy, grp = NULL, outcomes){
  x <- data.table(x)
  if (!comparator %in% unique(x[[strategy]])){
    stop("Chosen comparator strategy not in x")
  }
  if (is.null(grp)){
    grp = "grp"
    if (!"grp" %in% colnames(x)){
      x[, (grp) := 1] 
    }
  }
  indx_comparator <- which(x[[strategy]] == comparator)
  indx_treat <- which(x[[strategy]] != comparator)
  x_comparator <- x[indx_comparator]
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

# Survival quantiles -----------------------------------------------------------
#' Survival quantiles
#' 
#' Compute quantiles from survival curves.
#' @param x A \code{data.table} or \code{data.frame}.
#' @param probs A numeric vector of probabilities with values in \code{[0,1]}.
#' @param t A character scalar of the name of the time column.
#' @param surv_cols A character vector of the names of columns containing 
#' survival curves.
#' @param by A character vector of the names of columns to group by.
#' @return A \code{data.table} of quantiles of each survival curve in 
#' \code{surv_cols} by each group in \code{by}.
#' @examples 
#' library("data.table")
#' t <- seq(0, 10, by = .01)
#' surv1 <- seq(1, .3, length.out = length(t))
#' surv2 <- seq(1, .2, length.out = length(t))
#' strategies <- c("Strategy 1", "Strategy 2")
#' surv <- data.table(strategy = rep(strategies, each = length(t)),
#'                    t = rep(t, 2), 
#'                    surv = c(surv1, surv2))
#' surv_quantile(surv, probs = c(.4, .5), t = "t",
#'               surv_cols = "surv", by = "strategy")
#' @export
surv_quantile <- function (x, probs = .5, t, surv_cols, by) {
  x <- data.table(x)
  
  res <- vector(mode = "list", length = length(probs))
  for (i in 1:length(probs)){
    if (probs[i] > 1 | probs[i] < 0){
      stop("'prob' must be in the interval [0,1]",
           call. = FALSE)
    }  
    
    for (j in 1:length(surv_cols)){
      rows <- x[, .I[which(get(surv_cols[j]) <=  1 - probs[i])[1]], 
                by = by]$V1
      obs1_rows <- x[, .I[1], by = by]$V1
      na_rows <- which(is.na(rows))
      rows[is.na(rows)] <- obs1_rows[na_rows]
      quantile_name <- paste0("quantile_", surv_cols[j])
      if (j == 1){
        res[[i]] <- x[rows, c(by, t), with = FALSE]
        res[[i]][, ("prob") := probs[i]]
        setnames(res[[i]], t, quantile_name)
        setcolorder(res[[i]], c(by, "prob"))
      } else{
        res[[i]][, (quantile_name) := x[rows, t, with = FALSE]]
      }
      if (length(na_rows) > 0) {
        res[[i]][na_rows, (quantile_name) := NA]
      } 
    } # End loop over surv_cols
  } # End loop over probs
  res <- rbindlist(res)
  return(res)
}

# Summarize costs and QALYs ----------------------------------------------------
check_summarize <- function(x){
  if (is.null(x$costs_)) {
    stop("Cannot summarize costs without first simulating 'costs_' with '$sim_costs()'.",
         call. = FALSE)
  }
  
  if (is.null(x$qalys_)) {
    stop("Cannot summarize QALYs without first simulating 'qalys_' with '$sim_qalys()'.",
         call. = FALSE)
  }      
}

#' Summarize costs and effectiveness
#' 
#' Summarize costs and quality-adjusted life-years (QALYs) given output simulated
#' from an economic model. The summary output is used to perform 
#' cost-effectiveness analysis with [`cea()`] and [`cea_pw()`].
#' @param costs Simulated costs by category (objects of class [`costs`]). 
#' @param qalys Simulated QALYs (objects of class [`qalys`]).
#' @param by_grp If `TRUE`, then costs and QALYs are computed by subgroup. If
#' `FALSE`, then costs and QALYs are aggregated across all patients (and subgroups).
#' @details If mean costs and/or QALYs have already been computed 
#' (i.e., an average within a population), then there 
#' must be one observation for each discount rate (`dr`), 
#' PSA sample (`sample`), treatment strategy (`strategy_id`), 
#' and health state (`state_id`). Alternatively, there can be a column
#' denoting a patient (`patient_id`), in which case outcomes will first be
#' averaged across patients. A `grp_id` column can also be used so that
#' outcomes are computed for each subgroup (if `by_grp = TRUE`); otherwise it is assumed that 
#' there is only one subgroup.
#' @keywords internal
#' @return An object of class [`ce`].
summarize_ce <- function(costs, qalys, by_grp = FALSE) {
  patient_wt <- NULL
  by_cols <- c("dr", "sample", "strategy_id")
  if (by_grp) by_cols <- c(by_cols, "grp_id")
  
  summarize_wlos <- function(x, costs = TRUE, by_grp, by_cols){
    # Create grp ID column if missing
    if (by_grp == TRUE & !"grp_id" %in% colnames(x)){
      x[, ("grp_id") := 1]
    }
    
    # Some differences between cost and QALY output
    if (costs) {
      by_cols <- c("category", by_cols)
      sd_cols <- "costs"
    } else{
      sd_cols <- "qalys"
    }
    by_cols0 <- c(by_cols, "patient_id")
    if("patient_wt" %in% colnames(x)) by_cols0 <- c(by_cols0, "patient_wt")
    
    # Summarize
    if ("patient_id" %in% colnames(x)){ # Mean across patients
      x_summary <- x[, lapply(.SD, sum), by = by_cols0, .SDcols = sd_cols] 
      if ("patient_wt" %in% colnames(x)){ # Weighted mean
        x_summary <- x_summary[, lapply(.SD, stats::weighted.mean, w = patient_wt),
                               by = by_cols, .SDcols = sd_cols]
      } else{ # Non-weighted mean
        x_summary <- x_summary[, lapply(.SD, mean), by = by_cols, .SDcols = sd_cols]
      }
      
    } else{ # Mean already computed by health state, so sum across health states
      # Only for individual patient simulation
      x_summary <- x[, lapply(.SD, sum), by = by_cols, .SDcols = sd_cols]
    }
  }
  
  # QALYs
  qalys_summary <- summarize_wlos(qalys, costs = FALSE, by_grp, by_cols)
  
  # Costs
  costs_summary <- summarize_wlos(costs, costs = TRUE, by_grp, by_cols)
  costs_total <- costs_summary[, list(costs = sum(costs)), by = by_cols]
  costs_total[, ("category") := "total"]  
  costs_summary <- rbind(costs_summary, costs_total)
  
  # Combine
  ce <- list(costs = costs_summary, qalys = qalys_summary)
  if (by_grp == FALSE){
    ce[["costs"]][, ("grp_id") := 1]
    ce[["qalys"]][, ("grp_id") := 1]
  }
  class(ce) <- "ce"
  return(ce)
}

# Summary method ---------------------------------------------------------------
#' Summary method for cost-effectiveness object
#' 
#' Summarize a [`ce`] object by producing confidence intervals for quality-adjusted
#' life-years (QALYs) and each cost category with `summary.ce()` and format for 
#' pretty printing with `format.summary.ce()`.
#' @inheritParams icer
#' @param object A [`ce`] object. 
#' @param prob A numeric scalar in the interval `(0,1)` giving the confidence interval.
#' Default is 0.95 for a 95 percent interval. 
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' 
#' @return `summary.ce()` returns an object of class `summary.ce` that is a tidy
#'  `data.table` with the following columns:
#' 
#' \describe{
#' \item{dr}{The discount rate.}
#' \item{strategy}{The treatment strategy.}
#' \item{grp}{The patient subgroup.}
#' \item{type}{Either `"QALYs"` or `"Costs"`.}
#'  \item{category}{Category is always `"QALYs"` when `type == "QALYs"`; otherwise,
#'  it is the cost category.}
#' \item{estimate}{The point estimate computed as the average across the PSA samples.}
#' \item{lower}{The lower limit of the confidence interval.}
#' \item{upper}{The upper limit of the confidence interval.}
#' }
#' 
#' `format.summary.ce()` formats the table according to the arguments passed.
#' 
#' @details For an example, see [`IndivCtstm`].
#' @export
summary.ce <- function(object, prob = 0.95, labels = NULL,
                       ...) {
  qalys <- costs <- value <- dr <- grp <- order <- strategy <-  NULL
  alpha <- ci_alpha(prob)
  
  # Summarize the PSA separately for costs and QALYs
  res <- list(
    qalys = object$qalys[, list(type = "QALYs",
                                category = "QALYs",
                                estimate = mean(qalys),
                                lower = stats::quantile(qalys, alpha$lower),
                                upper = stats::quantile(qalys, alpha$upper)),
                         by = c("grp_id", "strategy_id", "dr")],
    costs = object$costs[, list(type = "Costs",
                                estimate = mean(costs),
                                lower = stats::quantile(costs, alpha$lower),
                                upper = stats::quantile(costs, alpha$upper)),
                         by = c("grp_id", "strategy_id", "dr", "category")]
  )

  # Combine costs and QALYs into a single table
  res <- rbindlist(res, use.names = TRUE)

  
  # Values of treatment strategies and subgroups
  set_labels(res, labels = labels)
  
  # Return
  setnames(res, c("grp_id", "strategy_id"), c("grp", "strategy"))
  setorderv(res, c("grp", "strategy", "dr"))
  setattr(res, "class", c("summary.ce", "data.table", "data.frame"))
  return(res[, ])
}

#' @rdname summary.ce
#' @inheritParams format.icer
#' @param x A `summary.ce` object.
#' @param dr_qalys Discount rate to subset to for quality-adjusted life-years (QALYs).
#' @param dr_costs Discount rate to subset to for costs.
#' @param pretty_names Logical. If `TRUE`, then the columns `strategy`, `grp`,
#'`outcome`, `dr`, and `value` are renamed (if they exist) to `Strategy`, `Group`,
#' `Outcome`, `Discount rate`, and `Value`.
#' @export
format.summary.ce <- function(x, digits_qalys = 2, digits_costs = 0,
                              dr_qalys = NULL, dr_costs = NULL,
                              pivot_from = "strategy", drop_grp = TRUE,
                              pretty_names = TRUE, ...) {
  dr <- type <- value <- estimate <- lower <- upper <- category <- 
    outcome <- grp <- NULL
  
  y <- copy(x)
  
  # Subset based on discount rate if specified
  if (!is.null(dr_qalys)) y <- y[dr %in% dr_qalys | type == "Costs"]
  if (!is.null(dr_costs)) y <- y[dr %in% dr_costs | type == "QALYs"]
    
  # Apply formatting
  y[, value := ifelse(type == "QALYs",
                      format_ci(estimate, lower, upper, costs = FALSE,
                                digits = digits_qalys),
                      format_ci(estimate, lower, upper, costs = TRUE,
                                digits = digits_costs))]
  y[, ("category") := ifelse(type == "QALYs",
                             category,
                             paste0(type, ": ", category))]
  setnames(y, "category", "outcome")
  y[, c("estimate", "lower", "upper", "type") := NULL]

  # Potentially pivot wider and drop groups
  id <- c("grp", "strategy", "dr", "outcome")
  y[, outcome := factor(outcome, levels = unique(y$outcome))]
  y <- format_summary_default(y, pivot_from = pivot_from,
                              id_cols = id,
                              drop_grp = drop_grp)
  
  # Pretty names
  if (pretty_names) {
    setnames(
      y, 
      c("grp", "strategy", "dr", "outcome", "value"),
      c("Group", "Strategy", "Discount rate", "Outcome", "Value"),
      skip_absent = TRUE
    )
  } 
  
  # Return
  return(y[, ])
}