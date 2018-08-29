context("cea.R unit tests")
rm(list = ls())

# Output for testing ----------------------------------------------------------
n_samples <- 1000

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c("Strategy 1, Grp 1", "Strategy 1, Grp 2", "Strategy 2, Grp 1",
              "Strategy 2, Grp 2", "Strategy 3, Grp 1", "Strategy 3, Grp 2")
c[[1]] <- rlnorm(n_samples, 2, .1)
c[[2]] <- rlnorm(n_samples, 2, .1)
c[[3]] <- rlnorm(n_samples, 11, .15)
c[[4]] <- rlnorm(n_samples, 11, .15)
c[[5]] <- rlnorm(n_samples, 11, .15)
c[[6]] <- rlnorm(n_samples, 11, .15)

# effectiveness
e <- c
e[[1]] <- rnorm(n_samples, 8, .2)
e[[2]] <- rnorm(n_samples, 8, .2)
e[[3]] <- rnorm(n_samples, 10, .8)
e[[4]] <- rnorm(n_samples, 10.5, .8)
e[[5]] <- rnorm(n_samples, 8.5, .6)
e[[6]] <- rnorm(n_samples, 11, .6)

# cost and effectiveness by strategy and simulation
ce <- data.table::data.table(sample = rep(seq(n_samples), length(e)),
                 strategy = rep(paste0("Strategy ", seq(1, 3)), 
                           each = n_samples * 2),
                 grp = rep(rep(c("Group 1", "Group 2"),
                               each = n_samples), 3),
                 cost = do.call("c", c), qalys = do.call("c", e))

# net benefits by willingess to pay
krange <- seq(100000, 120000, 500)
kval1 <- sample(krange, 1)
kval2 <- sample(krange, 1)
ce[, nb1 := qalys * kval1 - cost]
ce[, nb2 := qalys * kval2 - cost]

# Functions to use for testing ------------------------------------------------
# probabilistic sensitivity analysis
iceaR <- function(x, kval, grpname, output = c("mce", "evpi")){
  output <- match.arg(output)
  x <- x[grp == grpname] 
  x[, nb := qalys * kval - cost]
  nb <- data.table::dcast(x, sample ~ strategy, value.var = "nb")
  strategies <- seq(1, ncol(nb) - 1)
  nb_names <- paste0("nb", strategies)
  data.table::setnames(nb, colnames(nb), c("sample", nb_names))
  nb[, maxj := apply(nb[, -1, with = FALSE], 1, which.max)]
  nb[, maxj := factor(maxj, levels = strategies)]
  
  # mce
  if (output == "mce"){
    ret <- as.numeric(prop.table(table(nb$maxj)))
  }
  
  # evpi
  if (output == "evpi"){
    enb <- as.numeric(nb[, lapply(.SD, mean), .SDcols = nb_names])
    enb_maxj <- which.max(enb)
    nb[, nbpi := apply(nb[, 2:(ncol(nb) -1), with = FALSE], 1, max)]
    nb$nbci <- nb[[nb_names[enb_maxj]]]
    vpi <- nb$nbpi - nb$nbci
    ret <- mean(vpi)
  }
  
  # return
  return(ret)
}

# incemental changes
deltaR <- function(x, comparator, grpname){
  x <- x[grp == grpname]
  x.comparator <- x[strategy == comparator]
  x.treat <- x[strategy != comparator]
  x.treat.qalys <- data.table::dcast(x.treat, sample ~ strategy, value.var = c("qalys"))
  x.treat.cost <- data.table::dcast(x.treat, sample ~ strategy, value.var = c("cost"))
  delta.qalys <- x.treat.qalys[, -1, with = FALSE] - x.comparator$qalys
  delta.cost <- x.treat.cost[, -1, with = FALSE] - x.comparator$cost
  n_samples <- max(x$sample)
  ret <- data.table::data.table(sample = seq(1, n_samples), 
                                strategy = rep(unique(x.treat$strategy), each = n_samples),
                    grp = grpname,
                    ie = c(as.matrix(delta.qalys)), 
                    ic = c(as.matrix(delta.cost)))
  return(ret)
}

# ceac 
ceacR <- function(ix, kval, grpname) {
  ix <- ix[grp == grpname]
  ix[, nb := ie * kval - ic] 
  ceac <- ix[, .(prob = mean(nb >= 0)), by = "strategy"]
}

# Test icea function ----------------------------------------------------------
test_that("icea", {
  
  # function gets expected results
  icea_dt <-  icea(ce, k = krange, sample = "sample", strategy = "strategy",
                   grp = "grp", e = "qalys", c = "cost")
  kval <- sample(krange, 1)
  
  ## summary
  ### group 2
  ce_mean <- ce[grp == "Group 2", .(e_mean = mean(qalys), 
                                     c_mean = mean(cost)), by = "strategy"]
  ce_lower <- ce[grp == "Group 2", .(e_lower = quantile(qalys, .025), 
                                     c_lower = quantile(cost, .025)), by = "strategy"]
  ce_upper <- ce[grp == "Group 2", .(e_upper = quantile(qalys, .975), 
                                     c_upper = quantile(cost, .975)), by = "strategy"]
  summary_test <- data.table::data.table(strategy = ce_mean$strategy, 
                                        e_mean = ce_mean$e_mean,
                                        e_lower = ce_lower$e_lower, 
                                        e_upper = ce_upper$e_upper,
                                        c_mean = ce_mean$c_mean, 
                                        c_lower = ce_lower$c_lower,
                                        c_upper = ce_upper$c_upper)
  expect_equal(summary_test, icea_dt$summary[grp == "Group 2", -2, with = FALSE])
  
  # mce
  ### group 1
  mce <- icea_dt$mce[grp == "Group 1" &  k == kval]
  mce_test <- iceaR(ce, kval , "Group 1", output = "mce")
  expect_equal(mce$prob, mce_test)
  
  ### group 2
  mce <- icea_dt$mce[grp == "Group 2" &  k == kval]
  mce_test <- iceaR(ce, kval , "Group 2", output = "mce")
  expect_equal(mce$prob, mce_test)
  
  ## evpi
  ### group 1
  evpi <- icea_dt$evpi[grp == "Group 1" &  k == kval]
  evpi_test <- iceaR(ce, kval , "Group 1", output = "evpi")
  expect_equal(evpi$evpi, evpi_test)
  
  ## function works with other names
  ce2 = data.table::copy(ce)
  data.table::setnames(ce2, c("sample", "strategy", "grp"), c("samp", "strategy_name", "group"))
  icea_dt2 <-  icea(ce2, k = krange, sample = "samp", strategy = "strategy_name",
                   grp = "group", e = "qalys", c = "cost")
  evpi_v2 <- icea_dt2$evpi[group == "Group 1" &  k == kval]
  expect_equal(evpi_v2$evpi, evpi$evpi)
})

# Test icea_pw function -------------------------------------------------------
test_that("icea_pw", {
  
  ### function gets expected results
  icea_pw_dt <-  icea_pw(ce,  k = krange, comparator = "Strategy 1",
                         sample = "sample", strategy = "strategy", e = "qalys", c = "cost")
  kval <- sample(krange, 1)
  
  ## delta
  delta <- icea_pw_dt$delta
  delta_test <- deltaR(ce, comparator = "Strategy 1", grpname = "Group 1")
  expect_equal(delta[grp == "Group 1"], delta_test)
  
  ## summary
  # group 2
  delta_mean <- delta[grp == "Group 2", .(ie_mean = mean(ie), 
                                    ic_mean = mean(ic)), by = "strategy"]
  delta_lower <- delta[grp == "Group 2", .(ie_lower = quantile(ie, .025), 
                                     ic_lower = quantile(ic, .025)), by = "strategy"]
  delta_upper <- delta[grp == "Group 2", .(ie_upper = quantile(ie, .975), 
                                     ic_upper = quantile(ic, .975)), by = "strategy"]
  icer <- delta_mean$ic_mean/delta_mean$ie_mean
  summary_test <- data.table::data.table(strategy = delta_mean$strategy, 
                                         ie_mean = delta_mean$ie_mean,
                                         ie_lower = delta_lower$ie_lower, 
                                         ie_upper = delta_upper$ie_upper,
                                         ic_mean = delta_mean$ic_mean, 
                                         ic_lower = delta_lower$ic_lower,
                                         ic_upper = delta_upper$ic_upper, 
                                         icer = icer)
  expect_equal(summary_test, icea_pw_dt$summary[grp == "Group 2", -2, with = FALSE])
  
  ## ceac
  # group 1
  ceac <- icea_pw_dt$ceac[grp == "Group 1" & k == kval]
  ceac_test <- ceacR(delta, kval = kval, grpname = "Group 1")
  expect_equal(ceac$prob, ceac_test$prob)
  
  # group 2
  ceac <- icea_pw_dt$ceac[grp == "Group 2" & k == kval]
  ceac_test <- ceacR(delta, kval = kval, grpname = "Group 2")
  expect_equal(ceac$prob, ceac_test$prob)
  
  ## inmb
  # group 2
  inmb <- icea_pw_dt$inmb[k == kval & grp == "Group 2"]
  einmb_test <- delta[grp == "Group 2", .(einmb = mean(ie * kval - ic)), 
                     by = "strategy"]
  expect_equal(inmb$einmb, einmb_test$einmb)
  
  ### function works with other names
  ce2 = data.table::copy(ce)
  data.table::setnames(ce2, c("sample", "strategy", "grp"), c("samp", "strategy_name", "group"))
  icea_pw_dt2 <- icea_pw(ce2,  k = krange, comparator = "Strategy 1",
                         sample = "samp", strategy = "strategy_name", grp = "group",
                         e = "qalys", c = "cost")
  ceac_v2 <- icea_pw_dt2$ceac[group == "Group 2" & k == kval]
  expect_equal(ceac$prob, ceac_v2$prob)
})

# Test incr_effect function ---------------------------------------------------
test_that("incr_effect", {
  delta <- incr_effect(ce, comparator = "Strategy 1", sample = "sample", 
                       strategy = "strategy", grp = "grp",
                       outcomes = c("cost", "qalys"))
  expect_equal(delta[strategy == "Strategy 2" & grp == "Group 1", icost][1], 
               ce[strategy == "Strategy 2" & grp == "Group 1", cost][1] -
                 ce[strategy == "Strategy 1" & grp == "Group 1", cost][1])
  expect_equal(delta[strategy == "Strategy 2" & grp == "Group 2", icost][5], 
               ce[strategy == "Strategy 2" & grp == "Group 2", cost][5] -
                 ce[strategy == "Strategy 1" & grp == "Group 2", cost][5])
})