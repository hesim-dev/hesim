context("bootstrap.R unit tests")
library("flexsurv")

# Partitioned survival fits  ---------------------------------------------------
test_that("bootstrap.partsurvfit() works even some replications have errors ", {
  fit <- suppressWarnings(flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                     data = ovarian[1:2, ], dist="gamma"))
  fit <- partsurvfit(flexsurvreg_list(fit), data = ovarian[1:2, ])
  set.seed(101)
  expect_warning(bootstrap(fit, B = 5, max_errors = 10, silent = TRUE))
  set.seed(101)
})