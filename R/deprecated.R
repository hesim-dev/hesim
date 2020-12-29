# Deprecated functions ---------------------------------------------------------
#' Individualized cost-effectiveness analysis
#'
#' These functions are deprecated, use [cea()] and [cea_pw()] instead. 
#' @param x An object of simulation output characterizing the probability distribution
#' of clinical effectiveness and costs.?ic
#' @param ... Further arguments passed to or from other methods. 
#' @export
#' @rdname icea
icea <- function(x, ...) {
  .Deprecated("cea")
  UseMethod("cea")
}

#' @export
#' @rdname icea
icea_pw <- function(x, ...) {
  .Deprecated("cea_pw")
  UseMethod("cea_pw")
}

#' ICER table
#'
#' Generate a table of incremental cost-effectiveness ratios given output from 
#' [cea_pw()].
#'
#' @param x An object of class `cea_pw` returned by [cea_pw()].
#' @param k Willingness to pay.
#' @param cri If `TRUE`, credible intervals are computed; otherwise 
#' they are not.
#' @param prob A numeric scalar in the interval `(0,1)` giving the credible interval.
#' Default is 0.95 for a 95 percent credible interval. 
#' @param digits_qalys Number of digits to use to report QALYs.
#' @param digits_costs Number of digits to use to report costs.
#' @param output Should output be a `data.table` or a list of matrices for
#' each group.
#' @param rownames Row names for matrices when `output = "matrix"`.
#' @param colnames Column names for matrices when `output = "matrix"`.
#' @param drop If `TRUE`, then the result is coerced to the lowest possible dimension. 
#' Relevant if `output = "matrix"` and there is one group, in which case a single
#' matrix will be returned if `drop = TRUE` and a list of length 1 will be returned
#' if `drop = FALSE`.
#' @seealso [cea_pw()]
#' @return If `output = "matrix"`, then a list of matrices (or a matrix if
#' \code{drop = TRUE}) reporting incremental cost-effectiveness ratios (ICERs)
#' by group. Specifically, each matrix contains five rows for: (i) 
#' incremental quality-adjusted life-years (QALYs), (ii) incremental costs,
#' (iii) the incremental net monetary benefit (NMB), (iv) the ICER, 
#' and (v) a conclusion stating whether each strategy is cost-effective relative
#'  to a comparator. The number of columns is equal to the
#' number of strategies (including the comparator).
#' 
#' If `output = "data.table"`, then the results are reported as a `data.table`,
#' with one row for each strategy and group combination.
#' @export
icer_tbl <- function(x, k = 50000, cri = TRUE, prob = 0.95, 
                     digits_qalys = 2, 
                     digits_costs = 0, output = c("matrix", "data.table"),
                     rownames = NULL, colnames = NULL,
                     drop = TRUE){
  .Deprecated("icer")
  if (!inherits(x, "cea_pw")){
    stop("'x' must be an object of class 'cea_pw'",
         call. = FALSE)
  }
  if (prob > 1 | prob < 0){
    stop("'prob' must be in the interval (0,1)",
         call. = FALSE)
  }
  
  strategy <- attributes(x)$strategy
  grp <- attributes(x)$grp
  output <- match.arg(output)
  tbl <- copy(x$summary)
  tbl[, "icer" := get("ic_mean")/get("ie_mean")]  
  tbl[, "inmb_numeric" := k * get("ie_mean") - get("ic_mean")]
  
  # Formatting
  tbl[, "iqalys" := format_qalys(get("ie_mean"), digits = digits_qalys)]
  tbl[, "icosts" := format_costs(get("ic_mean"), digits = digits_costs)]
  tbl[, "icer" := format_costs(get("icer"), digits = digits_costs)]
  tbl[, "icer" := ifelse(get("ic_mean") < 0 & get("ie_mean") >= 0, "Dominates", get("icer"))]
  tbl[, "icer" := ifelse(get("ic_mean") > 0 & get("ie_mean") <= 0, "Dominated", get("icer"))]
  tbl[, "inmb" := format_costs(get("inmb_numeric"), digits = digits_costs)]
  
  if(cri){
    prob_lower <- (1 - prob)/2
    prob_upper <- 1 - prob_lower
    x$delta[, "inmb" := k * get("ie") - get("ic")]
    if (prob == 0.95){
      tbl[, "iqalys" := format_cri(get("iqalys"), get("ie_lower"), get("ie_upper"), 
                                   costs = FALSE,
                                   digits = digits_qalys)]
      tbl[, "icosts" := format_cri(get("icosts"), get("ic_lower"), get("ic_upper"),
                                   costs = TRUE,
                                   digits = digits_costs)]
      inmb_dt <- x$delta[, list(mean = mean(get("inmb")),
                                lower = stats::quantile(get("inmb"), prob_lower),
                                upper = stats::quantile(get("inmb"), prob_upper)),
                         by = c(strategy, grp)]
      tbl[, "inmb" := format_cri(get("inmb"), inmb_dt$lower, inmb_dt$upper,
                                 costs = TRUE,
                                 digits = digits_costs)]
    } else {
      cri_dt <- x$delta[, list(iqalys_lower = stats::quantile(get("ie"), prob_lower),
                               iqalys_upper = stats::quantile(get("ie"), prob_upper),
                               icosts_lower = stats::quantile(get("ic"), prob_lower),
                               icosts_upper = stats::quantile(get("ic"), prob_upper),
                               inmb_lower = stats::quantile(get("inmb"), prob_lower),
                               inmb_upper = stats::quantile(get("inmb"), prob_upper)),
                        by = c(strategy, grp)]
      tbl[, "iqalys" := format_cri(get("iqalys"), cri_dt$iqalys_lower, 
                                   cri_dt$iqalys_upper, costs = FALSE,
                                   digits = digits_qalys)]
      tbl[, "icosts" := format_cri(get("icosts"), cri_dt$icosts_lower, 
                                   cri_dt$icosts_upper, costs = TRUE,
                                   digits = digits_costs)]      
      tbl[, "inmb" := format_cri(get("inmb"), cri_dt$inmb_lower, 
                                 cri_dt$inmb_upper, costs = TRUE,
                                 digits = digits_costs)]       
    }
    x$delta[, "inmb" := NULL]
  } # end credible interval calculations
  tbl[, "conclusion" := ifelse(get("inmb_numeric") >= 0,
                               "Cost-effective", "Not cost-effective")]
  tbl <- tbl[, c(strategy, grp, "iqalys", "icosts", "inmb", "icer", "conclusion"),
             with = FALSE]
  
  if (output == "matrix"){
    tbl_list <- split(tbl, by = grp)
    mat_list <- vector(mode = "list", length = length(tbl_list))
    names(mat_list) <- names(tbl_list)
    n_strategies <- length(unique(tbl[[strategy]]))
    mat <- matrix(NA, nrow = 5, ncol = n_strategies + 1)
    if(is.null(rownames)){
      rownames(mat) <- c("Incremental QALYs", "Incremental costs", 
                         "Incremental NMB", "ICER", "Conclusion")
    } else{
      rownames(mat) <- rownames
    }
    comp_pos <- attributes(x)$comparator_pos
    if (is.null(colnames)){
      strategy_names <- rep(NA, ncol(mat))
      strategy_names[comp_pos] <- attributes(x)$comparator
      strategy_names[-comp_pos] <- as.character(tbl_list[[1]][[strategy]])
      colnames(mat) <- strategy_names
    } else{
      colnames(mat) <- colnames
    }
    for (i in 1:length(mat_list)){
      mat[1, -comp_pos] <- tbl_list[[i]]$iqalys
      mat[2, -comp_pos] <- tbl_list[[i]]$icosts
      mat[3, -comp_pos] <- tbl_list[[i]]$inmb
      mat[4, -comp_pos] <- tbl_list[[i]]$icer
      mat[5, -comp_pos] <- tbl_list[[i]]$conclusion
      mat[, comp_pos] <- "-"
      mat_list[[i]] <- mat 
    }
    if (drop){
      if(length(mat_list) == 1){
        mat_list <- mat_list[[1]]
      }
    }
    return(mat_list)
  } else{
    return(tbl)
  }
}

# Deprecated arguments ---------------------------------------------------------
deprecate_point_estimate <- function(old, new, is_new_missing) {
  if (!is.null(old)) { 
    warning("'point_estimate' is deprecated; use 'uncertainty' instead.",
            call. = FALSE)
  }
  if (!is.null(old) && (old == TRUE & is_new_missing == TRUE)) {
    return("none")
  } else{
    return(new)
  }
}

deprecate_bootstrap <- function(old, new, is_new_missing) {
  if (!is.null(old)) { 
    warning("'bootstrap' is deprecated; use 'uncertainty' instead.",
            call. = FALSE)
  }
  if (!is.null(old) && (old == TRUE & is_new_missing == TRUE)) {
    return("bootstrap")
  } else{
    return(new)
  }
}