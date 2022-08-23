context("test-ctstm.cpp unit tests")

tmat <- rbind(
  c(NA, 1, 2),
  c(NA, NA, 3),
  c(NA, NA, NA)
)

test_that("is_absorbing", {
  expect_equal(
    hesim:::C_test_is_absorbing(tmat),
    c(FALSE, FALSE, TRUE)
  )
})

test_that("trans_mat.trans_id()", {
  expect_equal(
    hesim:::C_test_trans_mat_trans_id(tmat, 0),
    c(0, 1)
  )
  expect_equal(
    hesim:::C_test_trans_mat_trans_id(tmat, 1),
    2
  )
  expect_equal(
    hesim:::C_test_trans_mat_trans_id(tmat, 2),
    integer()
  )
})

test_that("trans_mat.to()", {
  expect_equal(
    hesim:::C_test_trans_mat_to(tmat, 0),
    c(1, 2)
  )
  expect_equal(
    hesim:::C_test_trans_mat_to(tmat, 1),
    2
  )
  expect_equal(
    hesim:::C_test_trans_mat_to(tmat, 2),
    integer()
  )
})

test_that("trans_mat.n_trans_", {
  expect_equal(
    hesim:::C_test_trans_mat_n_trans(tmat),
    3
  )
})
