# Utility functions for plotting -----------------------------------------------
format_dollar <- function(x) {
  paste0("$", formatC(x, format = "d", big.mark = ","))
}

# Plot cost-effectiveness plane ------------------------------------------------
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
  if (!inherits(x, "cea_pw")){
    stop("'x' must be an object of class 'cea_pw'",
         call. = FALSE)
  }
  pdata <- copy(x$delta)
  
  strategy <- attr(x, "strategy")
  grp <- attr(x, "grp")
  
  # Passing custom names from user
  set_labels(pdata, labels = labels)
  
  # X and y limits
  xlim <- ceiling(max(x$delta[, ie]) * 1.2)
  ylim <- max(x$delta[, ic]) * 1.2
  
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