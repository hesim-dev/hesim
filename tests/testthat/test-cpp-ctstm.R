context("test-ctstm.cpp unit tests")

test_that("is_absorbing", {
  tmat <- rbind(c(NA, 1, 2), 
                c(NA, NA, 3), 
                c(NA, NA, NA))
  expect_equal(hesim:::C_ctstm_is_absorbing(tmat),
               c(FALSE, FALSE, TRUE))
})