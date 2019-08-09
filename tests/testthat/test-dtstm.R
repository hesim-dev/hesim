context("dtstm.R unit tests")
library("data.table")
rm(list = ls())

n_samples <- 3
strategies <- data.frame(strategy_id = c(1, 2))
n_strategies <- nrow(strategies)
patients <- data.frame(patient_id = 1)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)

trans1 <- c(200, 300, 500,
            0, 400, 100,
            0, 0, 600) 
trans1_mat <- matrix(trans1, nrow = 3, byrow = TRUE)
trans2 <- c(300, 300, 400,
            0, 550, 50,
            0, 0, 450)
trans2_mat <- matrix(trans2, nrow = 3, byrow = TRUE)
n_trans <- length(trans1)
n_states <- nrow(trans1_mat)

# transprob_tbl ----------------------------------------------------------------
test_that("transprob_tbl", {
  # Correct setup
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

# CohortDtstmTrans -------------------------------------------------------------
# Dirichlet distribution
transprobs <- transprob_tbl(data.frame(strategy_id = rep(1:2, each = n_trans),
                                       transition_id = rep(1:9, 2),
                                       alpha = c(trans1,trans2)),
                            dist = "dirichlet",
                            hesim_data = hesim_dat)
transmod_dirichlet <- create_CohortDtstmTrans(transprobs, n = n_samples)

## Tests
check_final_rows <- function(params){
  final_rows <- apply(params, 3, function(x) x[3, ])
  expect_true(all(final_rows[1, ] == 0))
  expect_true(all(final_rows[2, ] == 0))
  expect_true(all(final_rows[3, ] == 1))
}
check_final_rows(transmod_dirichlet$params)
expect_equal(dim(transmod_dirichlet$params)[3], 
             n_samples * n_strategies)

# "Custom" distribution
alpha <- rbind(trans1_mat, trans2_mat)
transprob_dist <- rdirichlet_mat(n = n_samples, alpha = alpha)
transprob_dist <- data.table(sample = rep(1:n_samples, 
                                          each = n_trans * 2),
                             patient_id = 1,
                             strategy_id = rep(rep(1:2, each = n_trans),
                                               n_samples),
                             transition_id = rep(1:n_trans, n_states * 2),
                             value = c(aperm(transprob_dist, 
                                             perm = c(2, 1, 3))))
transprobs <- transprob_tbl(transprob_dist,
                            dist = "custom")
transmod_custom <- create_CohortDtstmTrans(transprobs, n = n_samples)

## Tests
check_final_rows(transmod_custom$params)
for (v in c("sample", "strategy_id", "patient_id", "transition_id")){
  expect_true(all.equal(transmod_dirichlet[[v]], transmod_custom[[v]]))
}

# Fixed distribution
transprobs <- transprob_tbl(data.frame(strategy_id = rep(1:n_strategies, each = n_trans),
                                       transition_id = rep(1:9, 2),
                                       est = c(trans1,trans2)),
                            dist = "fixed",
                            hesim_data = hesim_dat)
transmod_fixed <- create_CohortDtstmTrans(transprobs, n = n_samples)

## Tests
expect_equal(dim(transmod_fixed$params)[3],
             n_strategies * n_samples)
expect_equal(transmod_fixed$params[,, 1], trans1_mat)
expect_equal(transmod_fixed$params[,, 4], trans2_mat)

