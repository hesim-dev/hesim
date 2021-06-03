context("params_mlogit.R unit tests")

# params_mlogit() works as expected --------------------------------------------
test_that("params_mlogit() works with a list of data frames", {
  p <- params_mlogit(
    coefs = list(
      sicker = data.frame(
        intercept = c(-0.33, -.2),
        treat = c(log(.75), log(.8))
      ),
      
      death = data.frame(
        intercept = c(-1, -1.2),
        treat = c(log(.6), log(.65))
      )
    )
  )
  expect_true(inherits(p, "params_mlogit"))
  expect_equal(p$n_samples, 2)
})

test_that("params_mlogit() works with a list of matrices", {
  p <- params_mlogit(
    coefs = list(
      sicker = matrix(1:4, 2, 2),
      death = matrix(1:4, 2, 2)
    )
  )
  expect_true(inherits(p, "params_mlogit"))
  expect_equal(colnames(p$coef[,, 1]), c("x1", "x2"))
})

test_that("params_mlogit() works with array inputs", {
  p <- params_mlogit(
    coefs = array(1:8, dim = c(2, 2, 2))
  )
  expect_true(inherits(p, "params_mlogit"))
  expect_equal(colnames(p$coef[,, 1]), c("x1", "x2"))
})

# params_mlogit() throws errors ------------------------------------------------
test_that("params_mlogit() requires numeric coefficients", {
  expect_error(
    params_mlogit(
      coefs = array(as.character(1:8), dim = c(2, 2, 2))
    ),
    "'coefs' must be a numeric 3D array."
  )
})

test_that("params_mlogit() coefficients must be a 3D array", {
  expect_error(
    params_mlogit(
      coefs = array(0, dim = rep(1, 4))
    ),
    "is_3d_array(x) is not TRUE",
    fixed = TRUE
  )
})

# summary.params_mlogit() ------------------------------------------------------
test_that("summary.params_mlogit()", {
  p <- params_mlogit(
    coefs = list(
      sicker = data.frame(
        intercept = c(-0.30, -.2),
        treat = c(-.3, -.1)
      ),
      
      death = data.frame(
        intercept = c(-1, -1.2),
        treat = c(-.5, -.4)
      )
    )
  )
  ps <- summary(p)
  expect_true(inherits(ps, "data.table"))
  expect_equal(ps$to, rep(c("sicker", "death"), each = 2))
  expect_equal(ps$term, rep(c("intercept", "treat"), 2))
  expect_equal(ps$estimate, c(-.25, -.2, -1.1, -.45))
})

# print.params_mlogit() ----------------------------------------------------------
test_that("print.params_mlogit() works as expected", {
  p <- params_mlogit(
    coefs = array(1:16, dim = c(4, 2, 2))
  )
  expect_output(print(p), "A \"params_mlogit\" object")
  expect_output(print(p), "Summary of coefficients:")
  expect_output(print(p), "Number of parameter samples: 4")
  expect_output(print(p), "Number of transitions: 2")
})
