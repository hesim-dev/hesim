context("params_mlogit_list.R unit tests")

p_ex <- params_mlogit_list(
  ## Transitions from sick state (sick -> sicker, sick -> death)
  sick = params_mlogit(
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
  ),
  
  ## Transitions from sicker state (sicker -> death)
  sicker = params_mlogit(
    coefs = list(
      death = data.frame(
        intercept = c(-1.5, -1.4),
        treat = c(log(.5), log(.55))
      )
    )
  )
)
# params_mlogit_list() works as expected ---------------------------------------
test_that("params_mlogit() works with a list of data frames", {
  expect_true(inherits(p_ex, "params_mlogit_list"))
  expect_equal(p_ex[[1]]$n_samples, 2)
})

# summary.params_mlogit_list() -------------------------------------------------
test_that("summary.params_mlogit_list()", {
  ps <- summary(p_ex)
  expect_true(inherits(ps, "data.table"))
  expect_equal(ps$from, c(rep("sick", 4), rep("sicker", 2)))
  expect_equal(ps$to, c(rep("sicker", 2), rep("death", 4)))
  expect_equal(ps$term, rep(c("intercept", "treat"), 3))
})

# print.params_mlogit_list() ---------------------------------------------------
test_that("print.params_mlogit_list()", {
  expect_output(print(p_ex), "A \"params_mlogit_list\" object")
  expect_output(print(p_ex), "Summary of coefficient estimates:")
  expect_output(print(p_ex), "Number of parameter samples: 2")
  expect_output(print(p_ex), "Number of transitions by starting state: 2 1")
})
