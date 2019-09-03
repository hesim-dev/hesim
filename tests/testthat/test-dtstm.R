context("dtstm.R unit tests")
library("data.table")
rm(list = ls())

# Helper functions -------------------------------------------------------------
sim_markov_chain <- function(p, x0, times, time_stop){
  n_states <- ncol(p[, , 1])
  n_times <- length(times)
  x <- matrix(NA, nrow = n_times, ncol = n_states)
  time_interval <- 1
  p_t <- p[, , 1]
  x[1, ] <- x0
  for (t in 2:n_times){
    if (times[t] > time_stop[time_interval]){
      time_interval <- time_interval + 1
      p_t <- p[, , time_interval]
    }
    x[t, ] <- x[t - 1, ] %*% p_t
  }
  rownames(x) <- times
  return(x)
}

as_trans_mats <- function(x){
  lapply(x, function (x) matrix(x, nrow = sqrt(length(x)),
                                byrow = TRUE))
}

to_probs <- function(x, vector = TRUE){
  res <- lapply(x, function (y) y/rowSums(y))
  if (vector){
    res <- unlist(lapply(res, t))
  }
  return(res)
}

n_states <- function(x) sqrt(length(x[[1]]))
n_trans <- function(x) length(x[[1]])

# test_sim_stateprobs(x, sample, strategy_id){
#   im <- x$input_mats
#   index <- which(im$strategy_id == strategy_id)
#   sample_env <- sample; strategy_id_env <- strategy_id
#   if (sample %in% colnames(x)){
#     x_i <- x[strategy_id == strategy_id_env] 
#     stprobs1 <- y$stateprobs_[strategy_id == strategy_id_env]
#   } else{
#     x_i <- x[sample == sample_env & strategy_id == strategy_id_env]
#     stprobs1 <- y$stateprobs_[sample == sample_env & strategy_id == strategy_id_env]
#   }
#   stprobs1 <- stprobs1[patient_id == 1]
#   stprobs1 <- matrix(stprobs1$prob, nrow = length(unique(stprobs1$t)))
#   p <- aperm(array(x_i$est, 
#                    dim = c(sqrt(length(x_i$est)),
#                            sqrt(length(x_i$est)),
#                            1)), # Adjust so 1 = length of times intervals
#              c(2, 1, 3))
#   x0 <- y$start_stateprobs
#   stprobs2 <- sim_markov_chain(p, x0 = y$start_stateprobs, 
#                                times = unique(y$stateprobs_$t), 
#                                time_stop = Inf) # Adjust time stop
#   expect_equal(stprobs1, stprobs2)
# }

# p <- transmod_dirichlet$params[,, 1, drop = FALSE]
# x0 <- c(1, rep(0, ncol(p) - 1))
# times <- seq(0, 2, 1)
# time_stop <- Inf
# sim_markov_chain(p, x0, times, time_stop)

# Data and parameters ----------------------------------------------------------
n_samples <- 3
strategies <- data.frame(strategy_id = c(1, 2))
n_strategies <- nrow(strategies)
patients <- data.frame(patient_id = 1)
n_patients <- nrow(patients)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)

# Transition probabilities vary by treatment strategy
alpha <- list(st1 = c(200, 300, 500,
                      0, 400, 100,
                      0, 0, 600),
              st2 = c(300, 300, 400,
                      0, 550, 50,
                      0, 0, 450))
alpha_mats <- as_trans_mats(alpha)

# transprob_tbl ----------------------------------------------------------------
test_that("transprob_tbl", {
  # Correct setup
  transprobs <- transprob_tbl(data.frame(strategy_id = rep(1:2, each = n_trans),
                                         transition_id = rep(1:9, 2),
                                         alpha = unlist(alpha)),
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

# CohortDtstmTrans (Dirichlet distribution) ------------------------------------
# Run
transprobs <- transprob_tbl(data.frame(strategy_id = rep(1:2, each = n_trans),
                                       transition_id = rep(1:9, 2),
                                       alpha = unlist(alpha)),
                            dist = "dirichlet",
                            hesim_data = hesim_dat)
transmod_dirichlet <- create_CohortDtstmTrans(transprobs, n = n_samples)
tmp <- transmod_dirichlet$sim_stateprobs(n_cycles = 2)$stateprobs_



## Test
check_final_rows <- function(params){
  final_rows <- apply(params, 3, function(x) x[3, ])
  expect_true(all(final_rows[1, ] == 0))
  expect_true(all(final_rows[2, ] == 0))
  expect_true(all(final_rows[3, ] == 1))
}
check_final_rows(transmod_dirichlet$params)
expect_equal(dim(transmod_dirichlet$params)[3], 
             n_samples * n_strategies * n_patients)

# CohortDtstmTrans ("Custom" distribution) ------------------------------------
# Run
transprob_dist <- rdirichlet_mat(n = n_samples, 
                                 alpha = do.call("rbind", alpha_mats))
transprob_dist <- data.table(sample = rep(1:n_samples, 
                                          each = n_trans(alpha) * 2),
                             patient_id = 1,
                             strategy_id = rep(rep(1:2, each = n_trans(alpha)),
                                               n_samples),
                             transition_id = rep(1:n_trans(alpha),
                                                 n_states(alpha) * 2),
                             value = c(aperm(transprob_dist, 
                                             perm = c(2, 1, 3))))
transprobs <- transprob_tbl(transprob_dist,
                            dist = "custom")
transmod_custom <- create_CohortDtstmTrans(transprobs, n = n_samples)

## Test
check_final_rows(transmod_custom$params)
for (v in c("sample", "strategy_id", "patient_id", "transition_id")){
  expect_true(all.equal(transmod_dirichlet[[v]], transmod_custom[[v]]))
}

# CohortDtstmTrans (Fixed distribution) ------------------------------------
# Run
transprobs <- transprob_tbl(data.frame(strategy_id = rep(1:n_strategies,
                                                         each = n_trans(alpha)),
                                       transition_id = rep(1:9, 2),
                                       est = to_probs(alpha_mats)),
                            dist = "fixed",
                            hesim_data = hesim_dat)
transmod_fixed <- create_CohortDtstmTrans(transprobs, n = n_samples)

## Test
expect_equal(dim(transmod_fixed$params)[3],
             n_strategies * n_samples * n_patients)
expect_equal(transmod_fixed$params[,, 1], to_probs(alpha_mats, FALSE)[[1]])
expect_equal(transmod_fixed$params[,, 4], to_probs(alpha_mats, FALSE)[[2]])

