context("tpmatrix unit tests")

# Transition probability matrix ------------------------------------------------
test_that("tpmatrix() works correctly", {
  p <- c(.7, .6)
  tpmat <- tpmatrix(
    C, p,
    0, 1
  )
  n <- length(p)
  expect_true(inherits(tpmat, "data.table"))
  expect_equal(tpmat$s1_s1, 1 - p)
  expect_equal(tpmat$s1_s2, p)
  expect_equal(tpmat$s2_s1, rep(0, n))
  expect_equal(tpmat$s2_s2, rep(1, n))
})

test_that("tpmatrix() works with complement argument", {
  pmat <- data.frame(s1_s1 = 0, s1_s2 = .5, s2_s1 = .3, s2_s2 = 0)

  # Character vector
  tpmat1 <- tpmatrix(pmat, complement = c("s1_s1", "s2_s2"))
  expect_equal(tpmat1$s1_s1, .5)
  expect_equal(tpmat1$s2_s2, .7)

  # Numeric vector
  tpmat2 <- tpmatrix(pmat, complement = c(1, 4))
  tpmat3 <- tpmatrix(pmat, complement = c(1L, 4L))
  expect_equal(tpmat1, tpmat2)
  expect_equal(tpmat1, tpmat3)
})

test_that("tpmatrix() works with states argument", {
  p <- tpmatrix(
    .5, .5, .5, .5,
    states = c("s1", "s2"), sep = "."
  )
  expect_equal(
    colnames(p),
    c("s1.s1", "s1.s2", "s2.s1", "s2.s2")
  )
})

test_that("tpmatrix() throws error if complement argument is incorrectly specified", {
  expect_error(
    tpmatrix(2, complement = data.frame(2)),
    "'complement' must either be a vector of integers or a character vector."
  )
})

test_that("tpmatrix() throws error if it is not a square msatrix", {
  expect_error(
    tpmatrix(1, 2, 3),
    "tpmatrix() must be a square matrix.",
    fixed = TRUE
  )
})

test_that("tpmatrix() throws error if states has wrong length", {
  expect_error(
    tpmatrix(1, 2, 3, 4, states = "s1"),
    paste0(
      "The length of 'states' must equal the square root of the number of ",
      "elements in the transition probability matrix."
    ),
    fixed = TRUE
  )
})

# Summarize transition probability matrix --------------------------------------
# Transition probability IDs
h <- hesim_data(
  strategies = data.table(strategy_id = 1:2),
  patients = data.table(patient_id = 1:3)
)
input_data <- expand(h, by = c("strategies", "patients"))
tpmat_id <- tpmatrix_id(input_data, n_samples = 2)

# Transition probability matrix
p_12 <- ifelse(tpmat_id$strategy_id == 1, .6, .7)
p <- tpmatrix(
  C, p_12,
  0, 1
)

test_that("summarize.tpmatrix() works without unflattening", {
  ps <- summary(p)
  expect_equal(colnames(ps), c("from", "to", "mean", "sd"))
  expect_equal(nrow(ps), 4)
  expect_equal(colnames(p), paste0(ps$from, "_", ps$to))
  expect_equivalent(ps$mean, apply(p, 2, mean))
})

test_that("summarize.tpmatrix() works with variables probs arguments", {
  ps <- summary(p, prob = .5)
  expect_equal(colnames(ps), c("from", "to", "mean", "sd", "50%"))

  ps <- summary(p, prob = c(.25, .75, .9))
  expect_equal(colnames(ps), c("from", "to", "mean", "sd", "25%", "75%", "90%"))
})

test_that("summarize.tpmatrix() works with unflattening", {
  ps <- summary(p, unflatten = TRUE)
  expect_equal(colnames(ps), c("mean", "sd"))
  expect_true(is.matrix(ps$mean[[1]]))
  states <- attr(p, "states")
  expect_equal(
    ps$mean[[1]],
    matrix(apply(p, 2, mean),
      nrow = 2, byrow = TRUE,
      dimnames = list(states, states)
    )
  )
})

test_that("summarize.tpmatrix() works with ID argument", {
  # Flattened
  ps <- summary(p, id = tpmat_id)
  expect_equal(nrow(input_data) * 4, nrow(ps))
  expect_equal(colnames(ps), c(
    "strategy_id", "patient_id",
    "from", "to", "mean", "sd"
  ))
  expect_true(all(ps$sd == 0))
  expect_true(all(ps[strategy_id == 1 & from == "s1" & to == "s2"]$mean == .6))
  expect_true(all(ps[strategy_id == 2 & from == "s1" & to == "s2"]$mean == .7))

  # Unflattened
  ps <- summary(p, id = tpmat_id, unflatten = TRUE)
  expect_equal(nrow(ps), nrow(input_data))
  expect_equal(
    unlist(ps[strategy_id == 1]$mean),
    rep(c(0.4, 0.0, 0.6, 1.0), times = 3)
  )
  expect_equal(
    unlist(ps[strategy_id == 2]$mean),
    rep(c(0.3, 0.0, 0.7, 1.0), times = 3)
  )
})

test_that("summarize.tpmatrix() works with ID argument and time intervals", {
  x <- expand(h, by = c("strategies", "patients"), times = c(0, 2))
  tpid <- tpmatrix_id(x, n_samples = 2)
  p2 <- tpmatrix(
    C, ifelse(tpid$strategy_id == 1, .6, .7),
    0, 1
  )

  # Flattened
  ps <- summary(p2, id = tpid, probs = .9)
  expect_equal(colnames(ps), c(
    "strategy_id", "patient_id",
    "time_id", "time_start", "time_stop",
    "from", "to", "mean", "sd", "90%"
  ))
  expect_equal(nrow(ps), nrow(x) * 4)

  # Unflattened
  ps <- summary(p2, id = tpid, unflatten = TRUE)
  expect_equal(nrow(ps), nrow(x))
  expect_true(all(unlist(ps$sd) == 0))
  expect_equal(
    unlist(ps[strategy_id == 2]$mean),
    rep(c(0.3, 0.0, 0.7, 1.0), times = 6)
  )
})

# Transition probability matrix IDs---------------------------------------------
strategies <- data.frame(strategy_id = c(1, 2))
patients <- data.frame(
  patient_id = seq(1, 3),
  patient_wt = c(1 / 2, 1 / 4, 1 / 4),
  gender = c("Female", "Female", "Male")
)
hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients
)

test_that("tpmatrix_id() returns errror if 'object' is not right class.", {
  expect_error(tpmatrix_id(2, 1))
})

test_that("tpmatrix_id() returns errror if 'object' is not expanded correctly.", {
  object <- expand(hesim_dat, by = c("strategies"))
  expect_error(
    tpmatrix_id(object, 1),
    paste0(
      "'object' must be expanded by 'strategy_id', 'patient_id',",
      " and optionally 'time_id'."
    )
  )
})

test_that("tpmatrix_id() returns correct numnber or rows and is the right class.", {
  input_data <- expand(hesim_dat, by = c("strategies", "patients"))
  tpmat_id <- tpmatrix_id(input_data, n_samples = 2)
  expect_equal(nrow(tpmat_id), nrow(input_data) * 2)
  expect_true(inherits(tpmat_id, "tpmatrix_id"))
})

test_that("tpmatrix_id() returns correct columns.", {
  input_data <- expand(hesim_dat,
    by = c("strategies", "patients"),
    times = c(0, 2)
  )
  tpmat_id <- tpmatrix_id(input_data, n_samples = 1)
  expect_equal(
    colnames(tpmat_id),
    c(
      "sample", "strategy_id", "patient_id", "patient_wt", "time_id",
      "time_start", "time_stop"
    )
  )
})

# Transition intensity matrix --------------------------------------------------
tmat <- rbind(
  c(NA, 1, 2),
  c(3, NA, 4),
  c(NA, NA, NA)
)
q12 <- c(.8, .7)
q13 <- c(.2, .3)
q21 <- c(.9, .8)
q23 <- c(1.1, 1.2)
q <- data.frame(q12, q13, q21, q23)
qmat <- qmatrix(q, trans_mat = tmat)

test_that("qmatrix() returns a 3D array", {
  expect_true(inherits(qmat, "array"))
  expect_equal(length(dim(qmat)), 3)
})

test_that("qmatrix() returns the correct diagonals", {
  expect_equal(mean(apply(qmat, 3, rowSums)), 0, tol = .00001)
})

# Transition intensity matrix (from msm) ---------------------------------------
set.seed(101)
library("msm")
qinit <- rbind(
  c(0, 0.28163, 0.01239),
  c(0, 0, 0.10204),
  c(0, 0, 0)
)
ptid <- sample(onc3p$patient_id, 200)
fit <- msm(state_id ~ time,
  subject = patient_id,
  data = onc3p[patient_id %in% ptid],
  covariates = ~ age + strategy_name, qmatrix = qinit
)

test_that("qmatrix.msm() works with factor covariates and 'newdata' is one row", {
  newdata <- data.frame(strategy_name = "New 1", age = 50)
  expect_equal(
    msm::qmatrix.msm(fit, newdata[1, , drop = FALSE], ci = "none"),
    qmatrix(fit, newdata, uncertainty = "none")[, , 1],
    check.attributes = FALSE
  )
})

test_that("qmatrix.msm() works with factor covariates and 'newdata' is multiple rows", {
  newdata <- data.frame(
    strategy_name = c("New 1", "New 2"),
    age = c(50, 55)
  )
  expect_equal(
    msm::qmatrix.msm(fit, newdata[2, , drop = FALSE], ci = "none"),
    qmatrix(fit, newdata, uncertainty = "none")[, , 2],
    check.attributes = FALSE
  )
})

test_that("qmatrix.msm() works with covariates that vary by transition", {
  fit <- update(fit, covariates = list("1-2" = ~ strategy_name + age))
  newdata <- data.frame(
    strategy_name = c("New 1", "New 2"),
    age = c(50, 55)
  )
  expect_equal(
    msm::qmatrix.msm(fit, newdata[2, , drop = FALSE], ci = "none"),
    qmatrix(fit, newdata, uncertainty = "none")[, , 2],
    check.attributes = FALSE
  )
})

test_that("qmatrix.msm() works with a hidden Markov model", {
  qinith <- rbind(
    c(0, exp(-6), exp(-9)),
    c(0, 0, exp(-6)),
    c(0, 0, 0)
  )
  hmod <- list(
    hmmNorm(mean = 100, sd = 16),
    hmmNorm(mean = 54, sd = 18),
    hmmIdent(999)
  )
  fith <- msm(fev ~ days,
    subject = ptnum,
    data = fev[fev$ptnum %in% 1:20, ],
    qmatrix = qinith,
    covariates = ~acute,
    hmodel = hmod,
    hcovariates = list(~acute, ~acute, NULL),
    hconstraint = list(acute = c(1, 1)),
    death = 3,
    method = "BFGS"
  )
  expect_equal(
    msm::qmatrix.msm(fith, covariates = list(acute = 0), ci = "none"),
    qmatrix(fith, data.frame(acute = 0), uncertainty = "none")[, , 1],
    check.attributes = FALSE
  )
})

test_that("qmatrix.msm() returns correct number of matrices with uncertainy = 'normal'", {
  newdata <- data.frame(strategy_name = c("New 1"), age = c(55))
  sim <- qmatrix(fit, newdata, uncertainty = "normal", n = 5)
  expect_true(dim(sim)[3] == 5)
})

test_that("qmatrix.msm() requires 'newdata' if covariates are included in the model.", {
  expect_error(
    qmatrix(fit),
    "'newdata' cannot be NULL if covariates are included in 'x'."
  )
})

test_that("qmatrix.msm() does not require 'newdata' if no covariates are included in the model.", {
  fit <- update(fit, covariates = ~1)
  expect_equal(
    msm::qmatrix.msm(fit, ci = "none"),
    qmatrix(fit, uncertainty = "none")[, , 1],
    check.attributes = FALSE
  )
})

# Matrix exponential -----------------------------------------------------------
test_that("expmat() returns an array where rows sum to 1.", {
  p <- expmat(qmat)
  expect_true(inherits(p, "array"))
  row_sums <- c(apply(p, 3, rowSums))
  expect_equal(mean(row_sums), 1, tolerance = .001, scale = 1)
})

test_that("expmat() works with matrix input.", {
  expect_true(inherits(expmat(qmat[, , 1]), "array"))
})

test_that("expmat() works with t as vector", {
  p <- expmat(qmat, t = c(1, 1))
  expect_equal(dim(p)[3], dim(qmat)[3] * 2)
  expect_equal(p[, , 1], p[, , 2])
  expect_equal(p[, , 3], p[, , 4])
})

test_that("expmat() is consisten with matrix multiplication ", {
  z <- diag(1, 3)
  p <- expmat(qmat, t = c(1, 2))
  expect_equal(
    z %*% p[, , 1] %*% p[, , 1],
    p[, , 2]
  )
})

test_that("expmat() returns error if x is not an array", {
  expect_error(
    expmat("Test error"),
    "'x' must be an array."
  )
})

# Convert to 3D array ----------------------------------------------------------
test_that("as_array3() returns a 3D array if 'x' is a square matrix", {
  expect_equal(
    dim(as_array3(matrix(1:16, 4, 4))),
    c(2, 2, 4)
  )
})

test_that("as_array3() throws error if 'x' is not a square matrix", {
  expect_error(
    as_array3(matrix(1:4, 2, 2)),
    "'x' must contain square matrices."
  )
})

# Relative risks ---------------------------------------------------------------
p_12 <- c(.7, .5)
p_23 <- c(.1, .2)
pmat <- as_array3(tpmatrix(
  C, p_12, .1,
  0, C, p_23,
  0, 0, 1
))
rr_12 <- runif(4, .8, 1)
rr_13 <- runif(4, .9, 1)
rr <- cbind(rr_12, rr_13)
pmat2 <- apply_rr(pmat, rr,
  index = list(c(1, 2), c(1, 3)),
  complement = list(c(1, 1), c(2, 2))
)

test_that("apply_rr() correctly multiplies relative risks", {
  expect_equal(pmat2[1, 2, ], rr_12 * pmat[1, 2, ])
  expect_equal(pmat2[1, 3, ], rr_13 * pmat[1, 3, ])
})

test_that("Row sums are correct with apply_rr()", {
  row_sums <- c(apply(pmat2, 3, rowSums))
  expect_equal(mean(row_sums), 1, tolerance = .001, scale = 1)
})

test_that("'index' argument in apply_rr() has correct dimensions", {
  expect_error(
    apply_rr(pmat, rr,
      index = list(c(1, 2), c(1, 3), c(2, 1)),
      complement = list(c(1, 1), c(2, 2))
    ),
    paste0(
      "'index' must contain the same number of matrix elements as the ",
      "number of columns in 'rr'."
    )
  )
})

test_that("'complement' argument in apply_rr() must have correct number of matrix elements", {
  expect_error(
    apply_rr(pmat, rr,
      index = list(c(1, 2), c(1, 3)),
      complement = list(c(1, 1), c(2, 2), c(3, 3), c(4, 4))
    ),
    paste0(
      "The number of matrix elements in 'complement' cannot be larger than the ",
      "number of rows in 'x'."
    )
  )
})

test_that("apply_rr() can only have one complementary column for each row in matrix", {
  expect_error(
    apply_rr(pmat, rr,
      index = list(c(1, 2), c(1, 3)),
      complement = list(c(1, 1), c(1, 2))
    ),
    "There can only be one complementary column in each row."
  )
})
