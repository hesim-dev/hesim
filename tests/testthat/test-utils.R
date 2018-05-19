context("utils.R unit tests")

# absorbing --------------------------------------------------------------------
test_that("absorbing", {
  tmat <- matrix(c(seq(1, 6), rep(NA, 3)), nrow = 3, ncol = 3, byrow = TRUE)
  expect_equal(absorbing(tmat), 3)
})

# list depth -------------------------------------------------------------------
list1 <- list(1)
list2 <- list(1, list1)
list3 <- list(1, list1, list2)

expect_equal(hesim:::list_depth(list1), 1)
expect_equal(hesim:::list_depth(list2), 2)
expect_equal(hesim:::list_depth(list3), 3)

# list to array ----------------------------------------------------------------
test_that("list_to_array", {
  # vector
  l1 <- list(a = c(1, 2, 3, 4), b = c(2, 3, 5, 10))
  a1 <- list_to_array(l1)
  expect_true(inherits(a1, "array"))
  expect_equal(dim(a1), c(1, 4, 2))
  expect_equal(l1$a, a1[,,1])
  
  # matrix
  l1 <- list(a = matrix(seq(1, 4), 2, 2),
             b = matrix(seq(5, 8), 2, 2))
  a1 <- list_to_array(l1)
  expect_true(inherits(a1, "array"))
  expect_equal(dim(a1), c(2, 2, 2))
  expect_equal(l1$a, a1[,,1])
  
  # non-matrix or non-vector
  obj <- list(a = data.frame(3))
  expect_error(list_to_array(obj))

})



