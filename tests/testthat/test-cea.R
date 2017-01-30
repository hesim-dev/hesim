context("Cost-effectiveness analysis")

# Output for testing ----------------------------------------------------------
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
ce <- data.table::data.table(sim = rep(seq(nsims), length(e)),
                 arm = rep(paste0("Arm ", seq(1, 3)), 
                           each = nsims * 2),
                 grp = rep(rep(c("Group 1", "Group 2"),
                               each = nsims), 3),
                 cost = do.call("c", c), qalys = do.call("c", e))

# net benefits by willingess to pay
krange <- seq(100000, 120000, 500)
kval1 <- sample(krange, 1)
kval2 <- sample(krange, 1)
ce[, nb1 := qalys * kval1 - cost]
ce[, nb2 := qalys * kval2 - cost]

# Functions to use for testing ------------------------------------------------
# probabilistic sensitivity analysis
psaR <- function(x, kval, grpname, output = c("mce", "evpi")){
  output <- match.arg(output)
  x <- x[grp == grpname] 
  x[, nb := qalys * kval - cost]
  nb <- data.table::dcast(x, sim ~ arm, value.var = "nb")
  arms <- seq(1, ncol(nb) - 1)
  nb.names <- paste0("nb", arms)
  data.table::setnames(nb, colnames(nb), c("sim", nb.names))
  nb[, maxj := apply(nb[, -1, with = FALSE], 1, which.max)]
  nb[, maxj := factor(maxj, levels = arms)]
  
  # mce
  if (output == "mce"){
    ret <- as.numeric(prop.table(table(nb$maxj)))
  }
  
  # evpi
  if (output == "evpi"){
    enb <- as.numeric(nb[, lapply(.SD, mean), .SDcols = nb.names])
    enb.maxj <- which.max(enb)
    nb[, nbpi := apply(nb[, 2:(ncol(nb) -1), with = FALSE], 1, max)]
    nb$nbci <- nb[[nb.names[enb.maxj]]]
    vpi <- nb$nbpi - nb$nbci
    ret <- mean(vpi)
  }
  
  # return
  return(ret)
}

# incemental changes
deltaR <- function(x, control, grpname){
  x <- x[grp == grpname]
  x.control <- x[arm == control]
  x.treat <- x[arm != control]
  x.treat.qalys <- data.table::dcast(x.treat, sim ~ arm, value.var = c("qalys"))
  x.treat.cost <- data.table::dcast(x.treat, sim ~ arm, value.var = c("cost"))
  delta.qalys <- x.treat.qalys[, -1, with = FALSE] - x.control$qalys
  delta.cost <- x.treat.cost[, -1, with = FALSE] - x.control$cost
  nsims <- max(x$sim)
  ret <- data.table::data.table(sim = seq(1, nsims), 
                    arm = rep(unique(x.treat$arm), each = nsims),
                    grp = grpname,
                    iqalys = c(as.matrix(delta.qalys)), 
                    icost = c(as.matrix(delta.cost)))
  return(ret)
}

# ceac 
ceacR <- function(ix, kval, grpname) {
  ix <- ix[grp == grpname]
  ix[, nb := iqalys * kval - icost] 
  ceac <- ix[, .(prob = mean(nb >= 0)), by = "arm"]
}

# Test psa function -----------------------------------------------------------
### defaults
psa.dt <-  psa(ce, k = krange, sim = "sim", arm = "arm",
                grp = "grp", e = "qalys", c = "cost")

test_that("psa", {
  kval <- sample(krange, 1)
  
  ## summary
  # group 2
  ce.mean <- ce[grp == "Group 2", .(qalys_mean = mean(qalys), 
                                     cost_mean = mean(cost)), by = "arm"]
  ce.lower <- ce[grp == "Group 2", .(qalys_lower = quantile(qalys, .025), 
                                     cost_lower = quantile(cost, .025)), by = "arm"]
  ce.upper <- ce[grp == "Group 2", .(qalys_upper = quantile(qalys, .975), 
                                     cost_upper = quantile(cost, .975)), by = "arm"]
  summary.test <- data.table::data.table(arm = ce.mean$arm, 
                             qalys_mean = ce.mean$qalys_mean,
                             qalys_lower = ce.lower$qalys_lower, 
                             qalys_upper = ce.upper$qalys_upper,
                              cost_mean = ce.mean$cost_mean, 
                             cost_lower = ce.lower$cost_lower,
                              cost_upper = ce.upper$cost)
  expect_equal(summary.test, psa.dt$summary[grp == "Group 2", -2, with = FALSE])
  
  ## mce
  # group 1
  mce <- psa.dt$mce[grp == "Group 1" &  k == kval]
  mce.test <- psaR(ce, kval , "Group 1", output = "mce")
  expect_equal(mce$prob, mce.test)
  
  # group 2
  mce <- psa.dt$mce[grp == "Group 2" &  k == kval]
  mce.test <- psaR(ce, kval , "Group 2", output = "mce")
  expect_equal(mce$prob, mce.test)
  
  ## evpi
  # group 1
  evpi <- psa.dt$evpi[grp == "Group 1" &  k == kval]
  evpi.test <- psaR(ce, kval , "Group 1", output = "evpi")
  expect_equal(evpi$evpi, evpi.test)
  
  # Group 2
  evpi <- psa.dt$evpi[grp == "Group 2" &  k == kval]
  evpi.test <- psaR(ce, kval , "Group 2", output = "evpi")
  expect_equal(evpi$evpi, evpi.test)
})

# Test psa_pw function --------------------------------------------------------
psa.pw.dt <-  psa_pw(ce,  k = krange, control = "Arm 1",
                     sim = "sim", arm = "arm", e = "qalys", c = "cost")

test_that("psa_pw", {
  kval <- sample(krange, 1)
  
  ## delta
  delta <- psa.pw.dt$delta
  delta.test <- deltaR(ce, control = "Arm 1", grpname = "Group 1")
  expect_equal(delta[grp == "Group 1"], delta.test)
  
  ## summary
  # group 2
  delta.mean <- delta[grp == "Group 2", .(iqalys_mean = mean(iqalys), 
                                    icost_mean = mean(icost)), by = "arm"]
  delta.lower <- delta[grp == "Group 2", .(iqalys_lower = quantile(iqalys, .025), 
                                     icost_lower = quantile(icost, .025)), by = "arm"]
  delta.upper <- delta[grp == "Group 2", .(iqalys_upper = quantile(iqalys, .975), 
                                     icost_upper = quantile(icost, .975)), by = "arm"]
  icer <- delta.mean$icost_mean/delta.mean$iqalys_mean
  summary.test <- data.table::data.table(arm = delta.mean$arm, 
                                         iqalys_mean = delta.mean$iqalys_mean,
                                         iqalys_lower = delta.lower$iqalys_lower, 
                                         iqalys_upper = delta.upper$iqalys_upper,
                                         icost_mean = delta.mean$icost_mean, 
                                         icost_lower = delta.lower$icost_lower,
                                         icost_upper = delta.upper$icost, 
                                         icer = icer)
  expect_equal(summary.test, psa.pw.dt$summary[grp == "Group 2", -2, with = FALSE])
  
  ## ceac
  # group 1
  ceac <- psa.pw.dt$ceac[grp == "Group 1" & k == kval]
  ceac.test <- ceacR(delta, kval = kval, grpname = "Group 1")
  expect_equal(ceac$prob, ceac.test$prob)
  
  # group 2
  ceac <- psa.pw.dt$ceac[grp == "Group 2" & k == kval]
  ceac.test <- ceacR(delta, kval = kval, grpname = "Group 2")
  expect_equal(ceac$prob, ceac.test$prob)
  
  ## einb
  # group 2
  einb <- psa.pw.dt$einb[k == kval & grp == "Group 2"]
  einb.test <- delta[grp == "Group 2", .(einb = mean(iqalys * kval - icost)), 
                     by = "arm"]
  expect_equal(einb$einb, einb.test$einb)

})