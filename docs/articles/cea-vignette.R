# Fitting a muti-state model with unknown transition times
## @knitr functions

## ---- DECISION ANALYSIS ------------------------------------------------------
## @knitr ce_output
nsims <- 1000

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c("Arm 1, Grp 1", "Arm 1, Grp 2", "Arm 2, Grp 1",
              "Arm 2, Grp 2", "Arm 3, Grp 1", "Arm 3, Grp 2")
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

# cost and effectiveness by arm and simulation
library("data.table")
ce <- data.table(sim = rep(seq(nsims), length(e)),
                             arm = rep(paste0("Arm ", seq(1, 3)), 
                                       each = nsims * 2),
                             grp = rep(rep(c("Group 1", "Group 2"),
                                           each = nsims), 3),
                             cost = do.call("c", c), qalys = do.call("c", e))
head(ce)

## @knitr enb_example
ce[, nb := 150000 * qalys - cost]
enb <- ce[, .(enb = mean(nb)), by = c("arm", "grp")]
enb <- dcast(enb, arm ~ grp, value.var = "enb")
print(enb)

## @knitr psa
library("hesim")
ktop <- 200000
psa.dt <-  psa(ce, k = seq(0, ktop, 500), sim = "sim", arm = "arm",
              grp = "grp", e = "qalys", c = "cost")

## @knitr psa_pw
psa.pw.dt <-  psa_pw(ce,  k = seq(0, ktop, 500), control = "Arm 1",
                     sim = "sim", arm = "arm", e = "qalys", c = "cost")

## @knitr mce_example_setup
set.seed(131)
library("knitr")
ce.nb <- dcast(ce[sim %in% sample(1:nsims, 10) & grp == "Group 2"], 
               sim ~ arm, value.var = "nb")
setnames(ce.nb, colnames(ce.nb), c("sim", "nb1", "nb2", "nb3"))
ce.nb[, maxj := apply(ce.nb[, .(nb1, nb2, nb3)], 1, which.max)]
ce.nb[, maxj := factor(maxj, levels = c(1, 2, 3))]

## @knitr mce_example
# library("DT")
# datatable(ce.nb, filter = "none", rownames = FALSE, 
#           options = list(searching = FALSE, scrollX = FALSE,
#                          bPaginate = FALSE, bInfo = FALSE)) %>% 
#             formatRound(colnames(ce.nb), 0)
kable(ce.nb, digits = 0, format = "html")
mce <- prop.table(table(ce.nb$maxj))
print(mce)

## @knitr mce_plot
library("ggplot2")
library("scales")
theme_set(theme_bw())
ggplot(psa.dt$mce, aes(x = k, y = prob, col = factor(arm))) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")

## @knitr evpi_example_a
armmax.g2 <- which.max(enb[[3]])
ce.nb[, nbpi := apply(ce.nb[, .(nb1, nb2, nb3)], 1, max)]
ce.nb[, nbci := ce.nb[[armmax.g2 + 1]]]
# datatable(ce.nb, filter = "none", rownames = FALSE,
#           options = list(searching = FALSE, scrollX = FALSE,
#                          bPaginate = FALSE, bInfo = FALSE)) %>% 
#   formatRound(colnames(ce.nb), 0)
kable(ce.nb, digits = 0, format = "html")

## @knitr evpi_example_b
enbpi <- mean(ce.nb$nbpi)
enbci <- mean(ce.nb$nbci)
print(enbpi)
print(enbci)
print(enbpi - enbci)

## @knitr evpi_plot
ggplot(psa.dt$evpi, aes(x = k, y = evpi)) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Expected value of perfect information") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")

## @knitr psa_custom
ce[, lys := qalys * 1.5]
cea.fun <- function(x) list(mean = mean(x), quant = quantile(x, c(.025, .975)))
psa.custom.dt <- psa(ce, k = seq(0, ktop, 500), sim = "sim", arm = "arm",
               grp = "grp", e = "qalys", c = "cost",
               custom_vars = c("cost", "lys", "qalys"), 
               custom_fun = cea.fun)

## @knitr outcome_dist1
psa.custom.dt$summary

## @knitr outcome_dist2
psa.custom.dt$custom.table

## @knitr ceplane_plot
head(psa.pw.dt$delta)
ylim <- max(psa.pw.dt$delta[, icost]) * 1.1
xlim <- ceiling(max(psa.pw.dt$delta[, iqalys]) * 1.1)
ggplot(psa.pw.dt$delta, aes(x = iqalys, y = icost, col = factor(arm))) + 
  geom_jitter(size = .5) + facet_wrap(~grp) + 
  xlab("Incremental QALYs") + ylab("Incremental cost") +
  scale_y_continuous(label = dollar, limits = c(-ylim, ylim)) +
  scale_x_continuous(limits = c(-xlim, xlim), breaks = seq(-6, 6, 2)) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm") +
  geom_abline(slope = 150000, linetype = "dashed") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

## @knitr ceac_plot
ggplot(psa.pw.dt$ceac, aes(x = k, y = prob, col = factor(arm))) +
  geom_line() + facet_wrap(~grp) + xlab("Willingess to pay") +
  ylab("Probability most cost-effective") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")

## @knitr icer
print(psa.pw.dt$summary)

## @knitr totevpi
w.dt <- data.table(grp = paste0("Group ", seq(1, 2)), w = c(0.25, .75))
evpi <- psa.dt$evpi
evpi <- merge(evpi, w.dt, by = "grp")
totevpi <- evpi[,lapply(.SD, weighted.mean, w = w),
                by = "k", .SDcols = c("evpi")]
ggplot(totevpi, aes(x = k, y = evpi)) +
  geom_line() + xlab("Willingess to pay") +
  ylab("Total EVPI") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), label = comma) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") + scale_colour_discrete(name = "Arm")

## @knitr totenb
ce <- merge(ce, w.dt, by = "grp")
totenb <- ce[, .(totenb = weighted.mean(nb, w = w)), by = c("arm")]

## @knitr evic1
ptenb.grp.max <- apply(as.matrix(enb[, -1]), 2, max)
ptenb.max <- sum(ptenb.grp.max * w.dt$w)
tenb.max <- max(totenb$totenb)
tnb <- c(ptenb.max, tenb.max)
names(tnb) <- c("Personalized total TENB", "One-size fits all TENB")

## @knitr evic2
evic <- tnb[1] - tnb[2]
names(evic) <- "EVIC"
print(evic)
print(evic/150000)
