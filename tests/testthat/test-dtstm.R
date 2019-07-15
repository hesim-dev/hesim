context("dtstm.R unit tests")
library("data.table")
rm(list = ls())

strategies <- data.frame(strategy_id = c(1, 2))
patients <- data.frame(patient_id = seq(1, 3))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)

# transprob_tbl ----------------------------------------------------------------
test_that("transprob_tbl", {
  # Correct setup
  trans1 <- c(200, 300, 500,
              0, 400, 100,
              0, 0, 600) 
  trans2 <- c(300, 300, 400,
              0, 550, 50,
              0, 0, 450)
  n_trans <- length(trans1)
  transprobs <- transprob_tbl(data.frame(strategy_id = rep(1:2, each = n_trans),
                                         transition_id = rep(1:9, 2),
                                         alpha = c(trans1,trans2)),
                              dist = "dirichlet",
                              hesim_data = hesim_dat)
  expect_equal(colnames(transprobs), 
               c("strategy_id", "transition_id", "alpha"))
  
  # Forget alpha column
  expect_error(transprob_tbl(data.frame(strategy_id = rep(1:2, each = n_trans),
                                         transition_id = rep(1:9, 2)),
                              dist = "dirichlet",
                              hesim_data = hesim_dat))
  
  # Incorrect number of rows
  expect_error(transprob_tbl(data.frame(strategy_id = c(rep(1, n_trans),
                                           rep(2, n_trans - 1)),
                                        transition_id = c(1:9, 2:9),
                                        alpha = c(trans1, trans2[-1])),
                              dist = "dirichlet",
                              hesim_data = hesim_dat))
  
  
  
})