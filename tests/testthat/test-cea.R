context("cea.R unit tests")
library("data.table")

# Output for testing ----------------------------------------------------------
n_samples <- 50

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c(
  "Strategy 1, Grp 1", "Strategy 1, Grp 2", "Strategy 2, Grp 1",
  "Strategy 2, Grp 2", "Strategy 3, Grp 1", "Strategy 3, Grp 2"
)
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
ce <- data.table(
  sample = rep(seq(n_samples), length(e)),
  strategy = rep(paste0("Strategy ", seq(1, 3)),
    each = n_samples * 2
  ),
  grp = rep(rep(c("Group 1", "Group 2"),
    each = n_samples
  ), 3),
  cost = do.call("c", c), qalys = do.call("c", e)
)

ce2 <- copy(ce)
setnames(ce2, c("sample", "strategy", "grp"), c("samp", "strategy_name", "group"))

# net benefits by willingess to pay
krange <- seq(0, 200000, by = 25000)
kval1 <- sample(krange, 1)
kval2 <- sample(krange, 1)
ce[, nmb1 := qalys * kval1 - cost]
ce[, nmb2 := qalys * kval2 - cost]

# Functions to use for testing ------------------------------------------------
# probabilistic sensitivity analysis
ceaR <- function(x, kval, grpname) {
  x <- x[grp == grpname]
  x[, nmb := qalys * kval - cost]
  nmb <- dcast(x, sample ~ strategy, value.var = "nmb")
  strategies <- seq(1, ncol(nmb) - 1)
  nmb_names <- paste0("nmb", strategies)
  setnames(nmb, colnames(nmb), c("sample", nmb_names))
  nmb[, maxj := apply(nmb[, -1, with = FALSE], 1, which.max)]
  nmb[, maxj := factor(maxj, levels = strategies)]
  prob_ce <- as.numeric(prop.table(table(nmb$maxj)))
  enmb <- as.numeric(nmb[, lapply(.SD, mean), .SDcols = nmb_names])
  enmb_maxj <- which.max(enmb)
  ret <- list()

  # mce
  ret$mce <- prob_ce

  # CEAF
  ret$ceaf <- prob_ce[enmb_maxj]

  # evpi
  nmb[, nmbpi := apply(nmb[, 2:(ncol(nmb) - 1), with = FALSE], 1, max)]
  nmb$nmbci <- nmb[[nmb_names[enmb_maxj]]]
  vpi <- nmb$nmbpi - nmb$nmbci
  ret$evpi <- mean(vpi)

  # return
  return(ret)
}

# incemental changes
deltaR <- function(x, comparator, grpname) {
  x <- x[grp == grpname]
  x.comparator <- x[strategy == comparator]
  x.treat <- x[strategy != comparator]
  x.treat.qalys <- dcast(x.treat, sample ~ strategy, value.var = c("qalys"))
  x.treat.cost <- dcast(x.treat, sample ~ strategy, value.var = c("cost"))
  delta.qalys <- x.treat.qalys[, -1, with = FALSE] - x.comparator$qalys
  delta.cost <- x.treat.cost[, -1, with = FALSE] - x.comparator$cost
  n_samples <- max(x$sample)
  ret <- data.table(
    sample = seq(1, n_samples),
    strategy = rep(unique(x.treat$strategy), each = n_samples),
    grp = grpname,
    ie = c(as.matrix(delta.qalys)),
    ic = c(as.matrix(delta.cost))
  )
  return(ret)
}

# ceac
ceacR <- function(ix, kval, grpname) {
  ix <- ix[grp == grpname]
  ix[, nmb := ie * kval - ic]
  ceac <- ix[, .(prob = mean(nmb >= 0)), by = "strategy"]
}

# Test cea() -------------------------------------------------------------------
cea_out <- cea(ce,
  k = krange, sample = "sample", strategy = "strategy",
  grp = "grp", e = "qalys", c = "cost"
)

test_that("cea() produced expected results", {

  # Function gets expected results
  kval <- sample(krange, 1)
  ceaR_1 <- ceaR(ce, kval, "Group 1")
  ceaR_2 <- ceaR(ce, kval, "Group 2")

  ## Summary
  ### Group 2
  ce_mean <- ce[grp == "Group 2", .(
    e_mean = mean(qalys),
    c_mean = mean(cost)
  ), by = "strategy"]
  ce_lower <- ce[grp == "Group 2", .(
    e_lower = quantile(qalys, .025),
    c_lower = quantile(cost, .025)
  ), by = "strategy"]
  ce_upper <- ce[grp == "Group 2", .(
    e_upper = quantile(qalys, .975),
    c_upper = quantile(cost, .975)
  ), by = "strategy"]
  summary_test <- data.table(
    strategy = ce_mean$strategy,
    e_mean = ce_mean$e_mean,
    e_lower = ce_lower$e_lower,
    e_upper = ce_upper$e_upper,
    c_mean = ce_mean$c_mean,
    c_lower = ce_lower$c_lower,
    c_upper = ce_upper$c_upper
  )
  expect_equal(summary_test, cea_out$summary[grp == "Group 2", -2, with = FALSE])

  # MCE
  ### Group 1
  mce <- cea_out$mce[grp == "Group 1" & k == kval]
  mce_test <- ceaR_1$mce
  expect_equal(mce$prob, mce_test)

  ### Group 2
  mce <- cea_out$mce[grp == "Group 2" & k == kval]
  mce_test <- ceaR_2$mce
  expect_equal(mce$prob, mce_test)

  ## CEAF
  ### Group 1
  ceaf <- cea_out$mce[best == 1 & grp == "Group 1" & k == kval]
  ceaf_test <- ceaR_1$ceaf
  expect_equal(ceaf$prob, ceaf_test)

  ### Group 2
  ceaf <- cea_out$mce[best == 1 & grp == "Group 2" & k == kval]
  ceaf_test <- ceaR_2$ceaf
  expect_equal(ceaf$prob, ceaf_test)

  ## EVPI
  ### Group 1
  evpi <- cea_out$evpi[grp == "Group 1" & k == kval]
  evpi_test <- ceaR_1$evpi
  expect_equal(evpi$evpi, evpi_test)

  # Function works with non-default names
  cea_out2 <- cea(ce2,
    k = krange, sample = "samp", strategy = "strategy_name",
    grp = "group", e = "qalys", c = "cost"
  )
  evpi_v2 <- cea_out2$evpi[group == "Group 1" & k == kval]
  expect_equal(evpi_v2$evpi, evpi$evpi)
})

# Test cea_pw() ----------------------------------------------------------------
cea_pw_out <- cea_pw(ce,
  k = krange, comparator = "Strategy 1",
  sample = "sample", strategy = "strategy", grp = "grp",
  e = "qalys", c = "cost"
)
cea_pw_out2 <- cea_pw(ce2,
  k = krange, comparator = "Strategy 1",
  sample = "samp", strategy = "strategy_name", grp = "group",
  e = "qalys", c = "cost"
)

test_that("cea_pw() produces expected results", {

  ### function gets expected results
  kval <- sample(krange, 1)

  ## delta
  delta <- cea_pw_out$delta
  delta_test <- deltaR(ce, comparator = "Strategy 1", grpname = "Group 1")
  expect_equal(delta[grp == "Group 1"], delta_test)

  ## summary
  # group 2
  delta_mean <- delta[grp == "Group 2", .(
    ie_mean = mean(ie),
    ic_mean = mean(ic)
  ), by = "strategy"]
  delta_lower <- delta[grp == "Group 2", .(
    ie_lower = quantile(ie, .025),
    ic_lower = quantile(ic, .025)
  ), by = "strategy"]
  delta_upper <- delta[grp == "Group 2", .(
    ie_upper = quantile(ie, .975),
    ic_upper = quantile(ic, .975)
  ), by = "strategy"]
  icer <- delta_mean$ic_mean / delta_mean$ie_mean
  summary_test <- data.table::data.table(
    strategy = delta_mean$strategy,
    ie_mean = delta_mean$ie_mean,
    ie_lower = delta_lower$ie_lower,
    ie_upper = delta_upper$ie_upper,
    ic_mean = delta_mean$ic_mean,
    ic_lower = delta_lower$ic_lower,
    ic_upper = delta_upper$ic_upper,
    icer = icer
  )
  expect_equal(summary_test, cea_pw_out$summary[grp == "Group 2", -2, with = FALSE])

  ## ceac
  # group 1
  ceac <- cea_pw_out$ceac[grp == "Group 1" & k == kval]
  ceac_test <- ceacR(delta, kval = kval, grpname = "Group 1")
  expect_equal(ceac$prob, ceac_test$prob)

  # group 2
  ceac <- cea_pw_out$ceac[grp == "Group 2" & k == kval]
  ceac_test <- ceacR(delta, kval = kval, grpname = "Group 2")
  expect_equal(ceac$prob, ceac_test$prob)

  ## inmb
  # group 2
  inmb <- cea_pw_out$inmb[k == kval & grp == "Group 2"]
  einmb_test <- delta[grp == "Group 2", .(einmb = mean(ie * kval - ic)),
    by = "strategy"
  ]
  expect_equal(inmb$einmb, einmb_test$einmb)

  # Function works with non-default names
  ceac_v2 <- cea_pw_out2$ceac[group == "Group 2" & k == kval]
  expect_equal(ceac$prob, ceac_v2$prob)
})

# Test icer_tbl() --------------------------------------------------------------
test_that("icer_tbl() produces expected results", {
  # No "credible interval"
  expect_warning(icer <- icer_tbl(cea_pw_out))
  expect_true(inherits(icer, "list"))
  expect_true(inherits(icer[[1]], "matrix"))

  # With "credible interval"
  expect_warning(icer <- icer_tbl(cea_pw_out, cri = TRUE))
  expect_true(inherits(icer, "list"))

  # data.table output
  expect_warning(icer <- icer_tbl(cea_pw_out, output = "data.table"))
  expect_true(inherits(icer, "data.table"))

  # Single group with output = "matrix"
  cea_pw_out2 <- cea_pw(ce[grp == "Group 1"],
    k = krange, comparator = "Strategy 1",
    sample = "sample", strategy = "strategy", e = "qalys", c = "cost"
  )
  expect_warning(icer <- icer_tbl(cea_pw_out2))
  expect_true(inherits(icer, "matrix"))
  expect_equal(ncol(icer), 3)

  # Single group with output = "matrix" and drop = FALSE
  expect_warning(icer <- icer_tbl(cea_pw_out2, drop = FALSE))
  expect_true(inherits(icer, "list"))

  # Specifying row and column names
  cols <- c("S1", "S2", "S3")
  rows <- c("iqalys", "icosts", "inmb", "icer", "conclusion")
  expect_warning(icer <- icer_tbl(cea_pw_out2, colnames = cols, rownames = rows))
  expect_true(inherits(icer, "matrix"))
  expect_equal(colnames(icer), cols)
  expect_equal(rownames(icer), rows)

  # Expect errors
  expect_error(expect_warning(icer_tbl(2)))
  expect_error(expect_warning(icer_tbl(cea_pw_out, prob = 1.4)))
})

# Test icer() ------------------------------------------------------------------
labs <- list(
  "strategy" = c(
    "s1" = "Strategy 1",
    "s2" = "Strategy 2",
    "s3" = "Strategy 3"
  ),
  "grp" = c(
    "g1" = "Group 1",
    "g2" = "Group 2"
  )
)

test_that("icer() and format.icer() return the correct columns", {
  # icer()
  x <- icer(cea_pw_out2)
  expect_equal(
    colnames(x),
    c("strategy", "grp", "outcome", "estimate", "lower", "upper")
  )

  # Formatting
  ## Pivoting
  ### Default
  y <- format(x)
  expect_true("Strategy 2" %in% colnames(y))

  ### No pivoting
  y <- format(x, pivot_from = NULL)
  expect_equal(colnames(y), c("Strategy", "Group", "Outcome", "Value"))

  ### Pivot strategy
  y <- format(x, pivot_from = "strategy")
  expect_true("Strategy 3" %in% colnames(y))

  ### Pivot group
  y <- format(x, pivot_from = "grp")
  expect_true("Group 1" %in% colnames(y))

  ### Pivot group and outcome
  y <- format(x, pivot_from = c("grp", "outcome"))
  expect_true(!"Outcome" %in% colnames(y))
  expect_true(!"Group" %in% colnames(y))
})

test_that("icer() correctly passes labels", {
  x <- icer(cea_pw_out, labels = labs)
  expect_equal(c("s2", "s3"), as.character(unique(x$strategy)))
  expect_equal(c("g1", "g2"), as.character(unique(x$grp)))
})

test_that("format.icer() will drop groups", {
  z <- cea_pw(ce[grp == "Group 1"],
    k = krange, comparator = "Strategy 1",
    sample = "sample", strategy = "strategy", grp = "grp",
    e = "qalys", c = "cost"
  )
  x <- format(icer(z))
  expect_true(!"grp" %in% colnames(x))
})

# Test plot_ceplane() ----------------------------------------------------------
p <- plot_ceplane(cea_pw_out, labels = labs)

test_that("plot_ceplane() returns ggplot", {
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_ceplane() correctly passes labels", {
  expect_equal(levels(p$data$strategy), names(labs$strategy))
  expect_equal(levels(p$data$grp), names(labs$grp))
})

test_that("plot_ceplane() throws error if 'x' is wrong class", {
  expect_error(
    plot_ceplane(2),
    "'x' must be of class 'cea_pw'."
  )
})

# Test plot_ceac() -------------------------------------------------------------
p1 <- plot_ceac(cea_out, labels = labs)
p2 <- plot_ceac(cea_pw_out, labels = labs)

test_that("plot_ceac returns ggplot from cea object", {
  expect_true(inherits(p1, "ggplot"))
})

test_that("plot_ceac returns ggplot from cea_pw object", {
  expect_true(inherits(p2, "ggplot"))
})

test_that("plot_ceac works with no labels", {
  expect_true(inherits(plot_ceac(cea_out), "ggplot"))
})

# Test plot_ceaf() -------------------------------------------------------------
p <- plot_ceaf(cea_out, labels = labs)

test_that("plot_ceaf() only uses data for optimal treatment strategy", {
  expect_equal(unique(p$data$best), 1)
})

test_that("plot_ceaf() throws error if 'x' is wrong class", {
  expect_error(
    plot_ceaf(2),
    "'x' must be of class 'cea'."
  )
})

# Test plot_evpi() -------------------------------------------------------------
test_that("plot_evpi() throws error if 'x' is wrong class", {
  expect_error(
    plot_evpi(2),
    "'x' must be of class 'cea'."
  )
})

test_that("plot_evpi() works with one group", {
  z <- cea(ce[grp == "Group 1"],
    k = krange, sample = "sample",
    strategy = "strategy", grp = "grp", e = "qalys", c = "cost"
  )
  expect_true(inherits(plot_evpi(z), "ggplot"))
})

# Test incr_effect() -----------------------------------------------------------
test_that("incr_effect", {
  # Default
  delta <- incr_effect(ce,
    comparator = "Strategy 1", sample = "sample",
    strategy = "strategy", grp = "grp",
    outcomes = c("cost", "qalys")
  )
  expect_equal(
    delta[strategy == "Strategy 2" & grp == "Group 1", icost][1],
    ce[strategy == "Strategy 2" & grp == "Group 1", cost][1] -
      ce[strategy == "Strategy 1" & grp == "Group 1", cost][1]
  )
  expect_equal(
    delta[strategy == "Strategy 2" & grp == "Group 2", icost][5],
    ce[strategy == "Strategy 2" & grp == "Group 2", cost][5] -
      ce[strategy == "Strategy 1" & grp == "Group 2", cost][5]
  )

  # Without group
  ce2 <- ce[grp == "Group 1"]
  ce2[, grp := NULL]
  delta <- incr_effect(ce2,
    comparator = "Strategy 1", sample = "sample",
    strategy = "strategy",
    outcomes = c("cost", "qalys")
  )
  expect_equal(
    delta[strategy == "Strategy 2", icost][3],
    ce[strategy == "Strategy 2", cost][3] -
      ce[strategy == "Strategy 1", cost][3]
  )
})
