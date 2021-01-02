# Utility functions for plotting -----------------------------------------------
format_dollar <- function(x) {
  paste0("$", formatC(x, format = "d", big.mark = ","))
}

# Cost-effectiveness plane -----------------------------------------------------
#' Plot cost-effectiveness plane
#' 
#' Plot a cost-effectiveness plane from the output of [`cea_pw()`] using [`ggplot2`]. 
#' Each point is a random draw of incremental costs (y-axis) and incremental QALYs (x-axis)
#' from a probabilistic sensitivity analysis.
#' @inheritParams set_labels
#' @param x A `cea_pw` object produced by [`cea_pw()`].
#' @param k Willingness to pay per QALY. 
#' @return A `ggplot` object.
#' @details See the [`cea_pw()`] documentation for an example. If there are multiple subgroups,
#' then a faceted plot is produced with one plot for each subgroup. 
#' @export 
plot_ceplane <- function(x, k = 50000, labels = NULL) {
  check_is_class(x, "cea_pw", "x")
  pdata <- copy(x$delta)
  
  # Some metadata
  strategy <- attr(x, "strategy")
  grp <- attr(x, "grp")
  
  # Passing custom names from user
  set_labels(pdata, labels = labels)
  
  # X and y limits
  xlim <- ceiling(max(x$delta[["ie"]]) * 1.2)
  ylim <- max(x$delta[["ic"]]) * 1.2
  
  # Main plot
  p <- ggplot2::ggplot(
    data = pdata,
    mapping = ggplot2::aes(x = .data[["ie"]], y = .data[["ic"]], 
                           col = .data[[strategy]])
  ) +
    ggplot2::geom_jitter(size = .5)  +
    ggplot2::xlab("Incremental QALYs") +
    ggplot2::ylab("Incremental costs") +
    ggplot2::scale_x_continuous(limits = c(-xlim, xlim),
                                breaks = pretty(seq(-xlim, xlim, length.out = 10), 
                                                n = 10)) +
    ggplot2::scale_y_continuous(limits = c(-ylim, ylim),
                                breaks = pretty(seq(-ylim, ylim, length.out = 10), 
                                                n = 10),
                                labels = format_dollar) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_colour_discrete(name = "Strategy") +
    ggplot2::geom_abline(slope = k, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0)
 
  # Add facets if more than one group
  n_grps <- length(unique(x$summary[[grp]]))
  if (n_grps > 1) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[grp]]))
  }
  
  # Return
  return(p)
}

# Cost-effectiveness acceptability curve ---------------------------------------
#' Plot cost-effectiveness acceptability curve
#' 
#' Plot a cost-effectiveness curve from either the output of [`cea()`] or
#' [`cea_pw()`] using [`ggplot2`]. The former compares all treatment strategies
#' simultaneously and uses the probabilistic sensitivity analysis (PSA) to compute
#' the probability that each strategy is the most cost-effective at a given 
#' willingness to pay value, while the latter uses the PSA to compute the probability
#' that each treatment is cost-effective relative to a comparator. 
#' 
#' @inheritParams set_labels
#' @param x An object of the appropriate class. 
#' @param ... Further arguments passed to and from methods. Currently unused. 
#' 
#' @details See the [`cea()`] documentation for an example. If there are multiple subgroups,
#' then a faceted plot is produced with one plot for each subgroup. 
#' @export
plot_ceac <- function(x, ...) {
  UseMethod("plot_ceac", x)
}

plot_ceac.default <- function(x, labels = NULL, ceaf = FALSE, ...) {
  best <- NULL
  
  # Get plotting data
  if (inherits(x, "cea_pw")){
    pdata <- copy(x$ceac)
  } else if (inherits(x, "cea")) {
    pdata <- copy(x$mce)
    if (ceaf) pdata <- pdata[best == 1]
  } else{
    x_class <- class(x)[1]
    stop(paste0("No default method for object of class ", x_class, "."))
  }
  
  # Some metadata
  strategy <- attr(x, "strategy")
  grp <- attr(x, "grp")
  max_k <- max(pdata$k)
  
  # Passing custom names from user
  set_labels(pdata, labels = labels)
  
  # Make plot
  p <- ggplot2::ggplot(
    data = pdata,
    mapping = ggplot2::aes(x = .data[["k"]], y = .data[["prob"]], 
                           col = .data[[strategy]])
  ) +
    ggplot2::facet_wrap(ggplot2::vars(.data[[grp]])) +
    ggplot2::geom_line() +
    ggplot2::xlab("Willingness to pay") +
    ggplot2::ylab("Probability most cost-effective") +
    ggplot2::scale_x_continuous(limits = c(0, max_k),
                                labels = format_dollar) + 
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, .2)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_colour_discrete(name = "Strategy")
  
  # Add facets if more than one group
  n_grps <- length(unique(x$summary[[grp]]))
  if (n_grps > 1) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[grp]]))
  }
  
  # Return
  return(p)
  
}

#' @export
#' @rdname plot_ceac
plot_ceac.cea_pw <- function(x, labels = NULL, ...) {
  plot_ceac.default(x, labels = labels)
}

#' @export
#' @rdname plot_ceac
plot_ceac.cea <- function(x, labels = NULL, ...) {
  plot_ceac.default(x, labels = labels)
}

# Cost-effectiveness acceptability frontier ------------------------------------
#' Plot cost-effectiveness acceptability frontier
#' 
#' Plot a cost-effectiveness acceptability frontier (CEAF) from the output of 
#' [`cea`] using [`ggplot2`]. The CEAF plots the probability
#' that the optimal treatment strategy (i.e., the strategy with the highest 
#' expected net monetary benefit) is cost-effective. 
#' @inheritParams set_labels
#' @param x A `cea` object produced by [`cea`].
#' @return A `ggplot` object.
#' @details See the [`cea()`] documentation for an example. If there are multiple subgroups,
#' then a faceted plot is produced with one plot for each subgroup. 
#' @export 
plot_ceaf <- function(x, labels = NULL) {
  check_is_class(x, "cea", "x")
  plot_ceac.default(x, labels = labels, ceaf = TRUE)
}

# Expected value of perfect information ----------------------------------------
#' Plot expected value of perfect information
#' 
#' Plot the expected value of perfect information (EVPI) from the output of 
#' [`cea()`] using [`ggplot2`]. Intuitively, the EVPI provides an estimate of the 
#' amount that a decision maker would be willing to pay to collect additional data 
#' and completely eliminate uncertainty.
#' @inheritParams set_labels
#' @param x A `cea` object produced by [`cea()`].
#' @return A `ggplot` object.
#' @details See the [`cea()`] documentation for an example. If there are multiple subgroups,
#' then a faceted plot is produced with one plot for each subgroup. 
#' @export 
plot_evpi <- function(x, labels = NULL) {
  check_is_class(x, "cea", "x")
  pdata <- copy(x$evpi)
  
  # Some metadata
  strategy <- attr(x, "strategy")
  grp <- attr(x, "grp")
  max_k <- max(pdata$k)
  
  # Passing custom names from user
  set_labels(pdata, labels = labels)
  
  # Main plot
  p <-  ggplot2::ggplot(
    data = pdata,
    mapping = ggplot2::aes(x = .data[["k"]], y = .data[["evpi"]])
  ) +
    ggplot2::geom_line() +
    ggplot2::xlab("Willingness to pay") +
    ggplot2::ylab("Expected value of perfect information") +
    ggplot2::scale_x_continuous(limits = c(0, max_k),
                                labels = format_dollar) +
    ggplot2::scale_y_continuous(breaks = pretty(pdata$evpi, n = 5),
                                labels = format_dollar) +
    ggplot2::theme(legend.position = "bottom")
  
  # Add facets if more than one group
  n_grps <- length(unique(x$summary[[grp]]))
  if (n_grps > 1) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[grp]]))
  }
  
  # Return
  return(p)
}
