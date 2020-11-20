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

test_that("qmatrix() returns a qmatrix object" , {
  expect_true(inherits(qmat, "qmatrix"))
})

test_that("qmatrix() returns the correct diagonals" , {
  expect_equal(-qmat[, 1], qmat[, 2] + qmat[, 3])
  expect_equal(-qmat[, 5], qmat[, 4] + qmat[, 6])
  expect_equal(-qmat[, 9], qmat[, 7] + qmat[, 9])
  expect_equal(c(0, 0), qmat[, 8])
})
