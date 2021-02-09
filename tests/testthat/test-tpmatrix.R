context("tpmatrix unit tests")

# Transition probability matrix ------------------------------------------------
test_that("tpmatrix() works correctly" , {
  p <- c(.7, .6)
  tpmat <- tpmatrix(
    C, p,
    0, 1
  )
  n <- length(p)
  expect_true(inherits(tpmat, "data.table"))
  expect_equal(tpmat$s1_s1, 1 - p)
  expect_equal(tpmat$s1_s2, p)
  expect_equal(tpmat$s2_s1,rep(0, n))
  expect_equal(tpmat$s2_s2,rep(1, n))
})

# Transition probability matrix IDs---------------------------------------------
strategies <- data.frame(strategy_id = c(1, 2))
patients <- data.frame(patient_id = seq(1, 3),
                       patient_wt = c(1/2, 1/4, 1/4),
                       gender = c("Female", "Female", "Male"))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)

test_that("tpmatrix_id() returns errror if 'object' is not right class." , {
  expect_error(tpmatrix_id(2, 1))
})

test_that("tpmatrix_id() returns errror if 'object' is not expanded correctly." , {
  object <- expand(hesim_dat, by = c("strategies"))
  expect_error(
    tpmatrix_id(object, 1),
    paste0("'object' must be either be expanded by 'strategy_id', 'patient_id',",
           " and optionally 'time_id'.")
  )
})

test_that("tpmatrix_id() returns correct numnber or rows and is the right class." , {
  input_data <- expand(hesim_dat, by = c("strategies", "patients"))
  tpmat_id <- tpmatrix_id(input_data, n_samples = 2)
  expect_equal(nrow(tpmat_id), nrow(input_data) * 2)
  expect_true(inherits(tpmat_id, "tpmatrix_id"))
})

test_that("tpmatrix_id() returns correct columns." , {
  input_data <- expand(hesim_dat, by = c("strategies", "patients"),
                       times = c(0, 2))
  tpmat_id <- tpmatrix_id(input_data, n_samples = 1)
  expect_equal(
    colnames(tpmat_id),
    c("sample", "strategy_id", "patient_id", "patient_wt", "time_id",
      "time_start", "time_stop")
  )             
})

# Transition intensity matrix --------------------------------------------------
tmat <- rbind(c(NA, 1, 2),
              c(3, NA, 4),
              c(NA, NA, NA)) 
q12 <- c(.8, .7)
q13 <- c(.2, .3)
q21 <- c(.9, .8)
q23 <- c(1.1, 1.2)
q <- data.frame(q12, q13, q21, q23)
qmat <- qmatrix(q, trans_mat = tmat)

test_that("qmatrix() returns a 3D array" , {
  expect_true(inherits(qmat, "array"))
  expect_equal(length(dim(qmat)), 3)
})

test_that("qmatrix() returns the correct diagonals" , {
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
fit <- msm(state_id ~ time, subject = patient_id, 
           data = onc3p[patient_id %in% ptid],
           covariates = ~ age + strategy_name, qmatrix = qinit)

test_that("qmatrix.msm() works with factor covariates and 'newdata' is one row" , {
  newdata <- data.frame(strategy_name = "New 1", age = 50)
  expect_equal(
    msm::qmatrix.msm(fit, newdata[1, , drop = FALSE], ci = "none"),
    qmatrix(fit, newdata, uncertainty = "none")[, , 1],
    check.attributes = FALSE
  )
})

test_that("qmatrix.msm() works with factor covariates and 'newdata' is multiple rows" , {
  newdata <- data.frame(strategy_name = c("New 1", "New 2"),
                        age = c(50, 55))
  expect_equal(
    msm::qmatrix.msm(fit, newdata[2, , drop = FALSE], ci = "none"),
    qmatrix(fit, newdata, uncertainty = "none")[, , 2],
    check.attributes = FALSE
  )
})

test_that("qmatrix.msm() works with covariates that vary by transition" , {
  fit <- update(fit, covariates = list("1-2" = ~ strategy_name + age))
  newdata <- data.frame(strategy_name = c("New 1", "New 2"),
                        age = c(50, 55))
  expect_equal(
    msm::qmatrix.msm(fit, newdata[2, , drop = FALSE], ci = "none"),
    qmatrix(fit, newdata, uncertainty = "none")[, , 2],
    check.attributes = FALSE
  )
})

test_that("qmatrix.msm() works with a hidden Markov model" , {
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
  fith <- msm(fev ~ days, subject = ptnum, 
              data = fev[fev$ptnum %in% 1:20, ],
              qmatrix = qinith, 
              covariates = ~acute,
              hmodel = hmod, 
              hcovariates = list(~ acute, ~ acute, NULL),
              hconstraint = list(acute = c(1, 1)), 
              death = 3,
              method = "BFGS")
  expect_equal(
    msm::qmatrix.msm(fith, covariates = list(acute = 0), ci = "none"),
    qmatrix(fith, data.frame(acute = 0), uncertainty = "none")[,,1],
    check.attributes = FALSE
  )
})  

test_that("qmatrix.msm() returns correct number of matrices with uncertainy = 'normal'" , {
  newdata <- data.frame(strategy_name = c("New 1"), age = c(55))
  sim <- qmatrix(fit, newdata, uncertainty = "normal", n = 5)
  expect_true(dim(sim)[3] == 5)
}) 

test_that("qmatrix.msm() requires 'newdata' if covariates are included in the model." , {
  expect_error(
    qmatrix(fit),
    "'newdata' cannot be NULL if covariates are included in 'x'."
  )
})

test_that("qmatrix.msm() does not require 'newdata' if no covariates are included in the model." , {
  fit <- update(fit, covariates =~ 1)
  expect_equal(
    msm::qmatrix.msm(fit, ci = "none"),
    qmatrix(fit, uncertainty = "none")[,,1],
    check.attributes = FALSE
  )
})
  
# Matrix exponential -----------------------------------------------------------
test_that("expmat() returns an array where rows sum to 1." , {
  p <- expmat(qmat)
  expect_true(inherits(p, "array"))
  row_sums <- c(apply(p, 3, rowSums))
  expect_equal(mean(row_sums), 1, tolerance = .001,scale = 1)
})

test_that("expmat() works with matrix input." , {
  expect_true(inherits(expmat(qmat[,,1]), "array"))
})

test_that("expmat() works with t as vector" , {
  p <- expmat(qmat, t = c(1, 1))
  expect_equal(dim(p)[3], dim(qmat)[3] * 2)        
  expect_equal(p[,, 1], p[,, 2])
  expect_equal(p[,, 3], p[,, 4])
})

test_that("expmat() is consisten with matrix multiplication " , {
  z <- diag(1, 3)
  p <- expmat(qmat, t = c(1, 2))
  expect_equal(
    z %*% p[,, 1] %*% p[, ,1],
    p[,, 2]
  )
})

test_that("expmat() returns error if x is not an array" , {
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
  0, C,     p_23,
  0, 0,     1
))
rr_12 <- runif(4, .8, 1)
rr_13 <- runif(4, .9, 1)
rr <- cbind(rr_12, rr_13)
pmat2 <- apply_rr(pmat, rr, 
                  index = list(c(1, 2), c(1, 3)),
                  complement = list(c(1, 1), c(2, 2)))

test_that("apply_rr() correctly multiplies relative risks" , {
  expect_equal(pmat2[1, 2, ], rr_12 * pmat[1, 2, ])
  expect_equal(pmat2[1, 3, ], rr_13 * pmat[1, 3, ])
})

test_that("Row sums are correct with apply_rr()" , {
  row_sums <- c(apply(pmat2, 3, rowSums))
  expect_equal(mean(row_sums), 1, tolerance = .001,scale = 1)
})

test_that("'index' argument in apply_rr() has correct dimensions" , {
  expect_error(
    apply_rr(pmat, rr, 
             index = list(c(1, 2), c(1, 3), c(2, 1)),
             complement = list(c(1, 1), c(2, 2))),
    paste0("'index' must contain the same number of matrix elements as the ",
           "number of columns in 'rr'.")
  )
})

test_that("'complement' argument in apply_rr() must have correct number of matrix elements" , {
  expect_error(
    apply_rr(pmat, rr, 
             index = list(c(1, 2), c(1, 3)),
             complement = list(c(1, 1), c(2, 2), c(3, 3), c(4, 4))),
    paste0("The number of matrix elements in 'complement' cannot be larger than the ",
            "number of rows in 'x'.")
  )
})

test_that("apply_rr() can only have one complementary column for each row in matrix" , {
  expect_error(
    apply_rr(pmat, rr, 
             index = list(c(1, 2), c(1, 3)),
             complement = list(c(1, 1), c(1,2))),
    "There can only be one complementary column in each row."
  )
})