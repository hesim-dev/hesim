# Fitting a muti-state model with unknown transition times
## @knitr functions

## ---- DECISION ANALYSIS ------------------------------------------------------
## @knitr mce_example_setup
library("data.table")
library("knitr")
out <- data.table(sim = seq(1:10),
                nb1 = rnorm(10, 10000, 1000), 
                nb2 = rnorm(10, 11000, 200),
                nb3 = rnorm(10, 9500, 2000))
out[, maxj := apply(out[, .(nb1, nb2, nb3)], 1, which.max)]
out[, maxj := factor(maxj, levels = c(1, 2, 3))]

## @knitr mce_example
kable(out, digits = 0)
mce <- prop.table(table(out$maxj))
print(mce)

## @knitr ce_output
n <- 1000

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c("Arm 1, Grp 1", "Arm 1, Grp 2", "Arm 2, Grp 1",
              "Arm 2, Grp 2", "Arm 3, Grp 1", "Arm 3, Grp 2")
c[[1]] <- rlnorm(n, 2, .1)
c[[2]] <- rlnorm(n, 2, .1)
c[[3]] <- rlnorm(n, 11, .15)
c[[4]] <- rlnorm(n, 11, .15)
c[[5]] <- rlnorm(n, 11, .15)
c[[6]] <- rlnorm(n, 11, .15)

# effectiveness
e <- c
e[[1]] <- rnorm(n, 8, .2)
e[[2]] <- rnorm(n, 8, .2)
e[[3]] <- rnorm(n, 10, .3)
e[[4]] <- rnorm(n, 10.5, .3)
e[[5]] <- rnorm(n, 9.5, .3)
e[[6]] <- rnorm(n, 11, .3)

# cost and effectiveness by arm and simulation
library("data.table")
ce <- data.table(sim = rep(seq(n), length(e)),
                             arm = rep(seq(1, 3), each = n * 2),
                             grp = rep(rep(c(1, 2), each = n), 3),
                             cost = do.call("c", c), qalys = do.call("c", e))

## @knitr enb_example
ce[, nb := 150000 * qalys - cost]
enb <- ce[, .(enb = mean(nb)), by = c("arm", "grp")]
enb <- dcast(enb, arm ~ grp, value.var = "enb")
print(enb)
ce[, nb := NULL]

## @knitr pcea
library("cea")
ce[, lys := qalys * 1.5]
cea.fun <- function(x) list(mean = mean(x), quant = quantile(x, c(.025, .975)))
ktop <- 200000
pcea.dt <-  pcea(ce, k = seq(0, ktop, 500), sim = "sim", arm = "arm",
                      grp = "grp", e = "qalys", c = "cost",
                      custom_vars = c("cost", "lys", "qalys"), custom_fun = cea.fun)

## @knitr mce_plot
library("ggplot2")
library("scales")
theme_set(theme_bw())
ggplot(pcea.dt$mce, aes(x = k, y = prob, col = factor(arm))) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")

## @knitr evpi_plot
ggplot(pcea.dt$evpi, aes(x = k, y = evpi)) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Expected value of perfect information") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")

## @knitr pcea_pw
pcea.pw.dt <-  pcea_pw(ce,  k = seq(0, ktop, 500), control = "1", e = "qalys", c = "cost",
                            custom_vars = c("cost", "lys", "qalys"))

## @knitr ceplane_plot
ylim <- max(pcea.pw.dt$delta[, icost]) * 2
xlim <- ceiling(max(pcea.pw.dt$delta[, iqalys]) * 1.5)
ggplot(pcea.pw.dt$delta, aes(x = iqalys, y = icost, col = factor(arm))) + 
  geom_jitter(size = .5) + facet_wrap(~grp) + 
  xlab("Incremental QALYs") + ylab("Incremental cost") +
  scale_y_continuous(label = dollar, limits = c(-ylim, ylim)) +
  scale_x_continuous(limits = c(-xlim, xlim), breaks = seq(-6, 6, 2)) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm") +
  geom_abline(slope = 150000, linetype = "dashed") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

## @knitr ceac_plot
ggplot(pcea.pw.dt$ceac, aes(x = k, y = prob, col = factor(arm))) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")

## @knitr icer
print(pcea.pw.dt$summary)

## @knitr totevpi
totevpi.dt <- totevpi(pcea.dt$evpi, grp = c(0, 1), w = c(0.25, .75))

