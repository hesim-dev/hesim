ggplot_tornado <- function(data, par = "par", par.values = "par.values",
                          labels = "labels", output.values = "output.values",
                          output.primary, output.label = "ICER"){
  p <- ggplot(data) +
    geom_segment(aes_string(y = par, yend = par,
                          x = output.values, xend = output.primary,
                     colour = labels), size = 3) +
    xlab(output.label) + ylab("") +
    geom_vline(xintercept = output.primary, linetype = 2, lwd = .1) +
    theme(legend.position = "bottom", legend.title = element_blank())
  return(p)
}

