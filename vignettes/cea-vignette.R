# Individualized cost-effectiveness analysis
## @knitr functions

## ---- DECISION ANALYSIS ------------------------------------------------------
## @knitr ce_output
nsims <- 1000

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c("Strategy 1, Grp 1", "Strategy 1, Grp 2", "Strategy 2, Grp 1",
              "Strategy 2, Grp 2", "Strategy 3, Grp 1", "Strategy 3, Grp 2")
c[[1]] <- rlnorm(nsims, 2, .1)
c[[2]] <- rlnorm(nsims, 2, .1)
c[[3]] <- rlnorm(nsims, 11, .15)
c[[4]] <- rlnorm(nsims, 11, .15)
c[[5]] <- rlnorm(nsims, 11, .15)
c[[6]] <- rlnorm(nsims, 11, .15)

# effectiveness
e <- c
e[[1]] <- rnorm(nsims, 8, .2)
e[[2]] <- rnorm(nsims, 8, .2)
e[[3]] <- rnorm(nsims, 10, .8)
e[[4]] <- rnorm(nsims, 10.5, .8)
e[[5]] <- rnorm(nsims, 8.5, .6)
e[[6]] <- rnorm(nsims, 11, .6)

# cost and effectiveness by strategy and simulation
library("data.table")
ce <- data.table(sim = rep(seq(nsims), length(e)),
                 strategy = rep(paste0("Strategy ", seq(1, 3)), 
                                each = nsims * 2),
                 grp = rep(rep(c("Group 1", "Group 2"),
                               each = nsims), 3),
                 cost = do.call("c", c), qalys = do.call("c", e))
head(ce)

## @knitr enmb_example
ce <- ce[, nmb := 150000 * qalys - cost]
enmb <- ce[, .(enmb = mean(nmb)), by = c("strategy", "grp")]
enmb <- dcast(enmb, strategy ~ grp, value.var = "enmb")
print(enmb)

## @knitr icea
library("hesim")
ktop <- 200000
icea.dt <-  icea(ce, k = seq(0, ktop, 500), sim = "sim", strategy = "strategy",
                 grp = "grp", e = "qalys", c = "cost")

## @knitr icea_pw
icea.pw.dt <-  icea_pw(ce,  k = seq(0, ktop, 500), comparator = "Strategy 1",
                       sim = "sim", strategy = "strategy", e = "qalys", c = "cost")

## @knitr mce_example_setup
set.seed(131)
library("knitr")
ce.nmb <- dcast(ce[sim %in% sample(1:nsims, 10) & grp == "Group 2"], 
                sim ~ strategy, value.var = "nmb")
setnames(ce.nmb, colnames(ce.nmb), c("sim", "nmb1", "nmb2", "nmb3"))
ce.nmb <- ce.nmb[, maxj := apply(ce.nmb[, .(nmb1, nmb2, nmb3)], 1, which.max)]
ce.nmb <- ce.nmb[, maxj := factor(maxj, levels = c(1, 2, 3))]

## @knitr mce_example
kable(ce.nmb, digits = 0, format = "html")
mce <- prop.table(table(ce.nmb$maxj))
print(mce)

## @knitr mce_plot
library("ggplot2")
library("scales")
theme_set(theme_bw())
ggplot(icea.dt$mce, aes(x = k, y = prob, col = factor(strategy))) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy")

## @knitr evpi_example_a
strategymax.g2 <- which.max(enmb[[3]])
ce.nmb <- ce.nmb[, nmbpi := apply(ce.nmb[, .(nmb1, nmb2, nmb3)], 1, max)]
ce.nmb <- ce.nmb[, nmbci := ce.nmb[[strategymax.g2 + 1]]]
kable(ce.nmb, digits = 0, format = "html")

## @knitr evpi_example_b
enmbpi <- mean(ce.nmb$nmbpi)
enmbci <- mean(ce.nmb$nmbci)
print(enmbpi)
print(enmbci)
print(enmbpi - enmbci)

## @knitr evpi_plot
ggplot(icea.dt$evpi, aes(x = k, y = evpi)) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Expected value of perfect information") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy")

## @knitr icea_summary
print(icea.dt$summary)

## @knitr icea_custom
ce <- ce[, lys := qalys * 1.5]
cea.fun <- function(x) list(mean = mean(x), quant = quantile(x, c(.025, .975)))
icea.custom.dt <- icea(ce, k = seq(0, ktop, 500), sim = "sim", strategy = "strategy",
                       grp = "grp", e = "qalys", c = "cost",
                       custom_vars = c("cost", "lys", "qalys"), 
                       custom_fun = cea.fun)

## @knitr outcome_dist1
icea.custom.dt$summary

## @knitr outcome_dist2
icea.custom.dt$custom.table

## @knitr ceplane_plot
head(icea.pw.dt$delta)
ylim <- max(icea.pw.dt$delta[, ic]) * 1.1
xlim <- ceiling(max(icea.pw.dt$delta[, ie]) * 1.1)
ggplot(icea.pw.dt$delta, aes(x = ie, y = ic, col = factor(strategy))) + 
  geom_jitter(size = .5) + facet_wrap(~grp) + 
  xlab("Incremental QALYs") + ylab("Incremental cost") +
  scale_y_continuous(label = dollar, limits = c(-ylim, ylim)) +
  scale_x_continuous(limits = c(-xlim, xlim), breaks = seq(-6, 6, 2)) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy") +
  geom_abline(slope = 150000, linetype = "dashed") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

## @knitr ceac_plot
ggplot(icea.pw.dt$ceac, aes(x = k, y = prob, col = factor(strategy))) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy")

## @knitr icer
print(icea.pw.dt$summary)

## @knitr totevpi
w.dt <- data.table(grp = paste0("Group ", seq(1, 2)), w = c(0.25, .75))
evpi <- icea.dt$evpi
evpi <- merge(evpi, w.dt, by = "grp")
totevpi <- evpi[,lapply(.SD, weighted.mean, w = w),
                by = "k", .SDcols = c("evpi")]
ggplot(totevpi, aes(x = k, y = evpi)) +
  geom_line() + xlab("Willingess to pay") +
  ylab("Total EVPI") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Strategy")

## @knitr totenmb
ce <- merge(ce, w.dt, by = "grp")
totenmb <- ce[, .(totenmb = weighted.mean(nmb, w = w)), by = c("strategy")]

## @knitr evic1
ptenmb.grp.max <- apply(as.matrix(enmb[, -1]), 2, max)
ptenmb.max <- sum(ptenmb.grp.max * w.dt$w)
tenmb.max <- max(totenmb$totenmb)
tnmb <- c(ptenmb.max, tenmb.max)
names(tnmb) <- c("Personalized total TENMB", "One-size fits all TENMB")

## @knitr evic2
evic <- tnmb[1] - tnmb[2]
names(evic) <- "EVIC"
print(evic)
print(evic/150000)