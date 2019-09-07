context("dtstm.R unit tests")
library("data.table")
rm(list = ls())

# Helper functions -------------------------------------------------------------
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

test_sim_stateprobs <- function(x, sample_val = 1, strategy_id_val = 1, 
                                patient_id_val = 1){
  # Simulations with hesim
  stprobs1 <- x$stateprobs_[sample == sample_val &
                              strategy_id == strategy_id_val &
                              patient_id == patient_id_val]
  stprobs1 <-  matrix(stprobs1$prob, nrow = length(unique(stprobs1$t)))
  
  # Simulate Markov chain with R function
  index <- which(x$params$sample == sample_val &
                   x$params$strategy_id == strategy_id_val &
                   x$params$patient_id == patient_id_val)
  p <- x$params$value[,, index, drop = FALSE]
  n_states <- ncol(p[,, 1])
  time_stop <- x$params$time_intervals$time_stop
  if (is.null(time_stop)) time_stop <- Inf
  stprobs2 <- sim_markov_chain(p = p,
                               x0 = c(1, rep(0, n_states - 1)),
                               times = unique(x$stateprobs_$t),
                               time_stop = time_stop)
  # Test
  expect_equal(c(stprobs1), c(stprobs2))
}

# Final rows are [0 0 ... 1]
check_final_rows <- function(p_array){
  final_rows <- apply(p_array, 3, function(x) x[3, ])
  n_states <- ncol(p_array[,, 1])
  for (i in 1:(n_states - 1)){
    expect_true(all(final_rows[i, ] == 0))
  }
  expect_true(all(final_rows[n_states, ] == 1))
}

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
transmod_dirichlet$sim_stateprobs(n_cycles = 2)

## Test
### Final rows are [0 0 ... 1]
check_final_rows(transmod_dirichlet$params$value)
expect_equal(dim(transmod_dirichlet$params$value)[3],
             n_samples * n_strategies * n_patients)

### Simulated state probabilities are correct
test_sim_stateprobs(transmod_dirichlet)
test_sim_stateprobs(transmod_dirichlet, sample_val = 2)
test_sim_stateprobs(transmod_dirichlet, sample_val = 3, strategy_id_val = 2)

# CohortDtstmTrans ("Custom" distribution) -------------------------------------
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
transmod_custom$sim_stateprobs(n_cycles = 4)

## Test
### Final rows are [0 0 ... 1]
check_final_rows(transmod_custom$params$value)

### ID values are correct
for (v in c("sample", "strategy_id", "patient_id")){
  expect_true(all.equal(transmod_dirichlet$params[[v]], 
                        transmod_custom$params[[v]]))
}

### Simulated state probabilities are correct
test_sim_stateprobs(transmod_custom, sample_val = 2, strategy_id_val = 2)
test_sim_stateprobs(transmod_custom, sample_val = 3, strategy_id_val = 1)

# CohortDtstmTrans (Fixed distribution) ------------------------------------
# Run
transprobs <- transprob_tbl(data.frame(strategy_id = rep(1:n_strategies,
                                                         each = n_trans(alpha)),
                                       transition_id = rep(1:9, 2),
                                       est = to_probs(alpha_mats)),
                            dist = "fixed",
                            hesim_data = hesim_dat)
transmod_fixed <- create_CohortDtstmTrans(transprobs, n = n_samples)
transmod_fixed$sim_stateprobs(n_cycles = 3)

## Test
### Transition probability array is correct
expect_equal(dim(transmod_fixed$params$value)[3],
             n_strategies * n_samples * n_patients)
expect_equal(transmod_fixed$params$value[,, 1], to_probs(alpha_mats, FALSE)[[1]])
expect_equal(transmod_fixed$params$value[,, 4], to_probs(alpha_mats, FALSE)[[2]])

### Final rows are [0 0 ... 1]
check_final_rows(transmod_fixed$params$value)

### Simulated state probabilities are correct
test_sim_stateprobs(transmod_fixed)
test_sim_stateprobs(transmod_fixed, sample_val = 2)
test_sim_stateprobs(transmod_fixed, sample_val = 3, strategy_id_val = 2)

