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

# Matrix exponential -----------------------------------------------------------
test_that("expmat() returns an array where rows sum to 1." , {
  p <- expmat(qmat)
  expect_true(inherits(p, "array"))
  row_sums <- c(apply(p, 3, rowSums))
  expect_equal(mean(row_sums), 1, tolerance = .001,scale = 1)
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