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

test_sim_stateprobs <- function(x, sample_val = 1, strategy_id_val = 1, 
                                patient_id_val = 1){
  # Simulations with hesim
  stprobs1 <- x$stateprobs_[sample == sample_val &
                              strategy_id == strategy_id_val &
                              patient_id == patient_id_val]
  stprobs1 <-  matrix(stprobs1$prob, nrow = length(unique(stprobs1$t)))
  
  # Simulate Markov chain with R function
  tm <- x$trans_model
  index <- which(tm$params$sample == sample_val &
                 tm$params$strategy_id == strategy_id_val &
                 tm$params$patient_id == patient_id_val)
  p <- tm$params$value[,, index, drop = FALSE]
  n_states <- ncol(p[,, 1])
  time_stop <- tm$params$time_intervals$time_stop
  if (is.null(time_stop)) time_stop <- Inf
  stprobs2 <- sim_markov_chain(p = p,
                               x0 = c(1, rep(0, n_states - 1)),
                               times = unique(x$stateprobs_$t),
                               time_stop = time_stop)
  # Test
  expect_equal(c(stprobs1), c(stprobs2))
}

apply_rr <- function(x, rr){
  x[upper.tri(x)] <- x[upper.tri(x)] * rr
  for (i in 1:(nrow(x) - 1)){
    x[i, i] <- 1 - sum(x[i, (i + 1):ncol(x)])
  }
  return(x)
}

make_alpha <- function(x, rr = 1, n = c(500, 800, 700)){
  return(apply_rr(x, rr) * n)
}

# Data and parameters ----------------------------------------------------------
n_samples <- 3
strategies <- data.frame(strategy_id = c(1, 2))
n_strategies <- nrow(strategies)
patients <- data.frame(patient_id = 1:2)
n_patients <- nrow(patients)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)
time_start <- c(0, 5)
n_times <- length(time_start)
n_states <- 3
tprob <- matrix(c(.2, .3, .5,
                   0, .8, .2,
                   0, 0, 1),
                ncol = 3, nrow = 3, byrow = TRUE)
tprob_array <- array(NA, dim = c(n_states, n_states, 2, n_patients,
                                 n_strategies, n_samples))

# Time ID = 1
tprob_array[, , 1, 1, 1, ] <- rdirichlet_mat(n_samples, make_alpha(tprob))
tprob_array[, , 1, 1, 2, ] <- rdirichlet_mat(n_samples, make_alpha(tprob, .7))
tprob_array[, , 1, 2, 1, ] <- rdirichlet_mat(n_samples, make_alpha(tprob, .9))
tprob_array[, , 1, 2, 2, ] <- rdirichlet_mat(n_samples, make_alpha(tprob, .9))

# Time ID = 2
tprob_array[, , 2, 1, 1, ] <- tprob_array[, , 1, 1, 1, ] 
tprob_array[, , 2, 1, 2, ] <- tprob_array[, , 1, 1, 1, ] 
tprob_array[, , 2, 2, 1, ] <- tprob_array[, , 1, 2, 1, ]
tprob_array[, , 2, 2, 2, ] <- tprob_array[, , 1, 2, 1, ]

tprob_array <- aperm(tprob_array, perm = c(6:3, 1, 2))

# tparams_transprobs.array -----------------------------------------------------
params_tprob <- tparams_transprobs(tprob_array, times = time_start)
params_tprob_t1 <- tparams_transprobs(tprob_array[,,,1,,, drop = FALSE])

test_that("Extra arguments with tparams_transprobs.array" , {
  expect_equal(tparams_transprobs(tprob_array, times = time_start, grp_id = 1)$grp_id,
               rep(1, prod(dim(tprob_array)[1:4])))
  expect_equal(tparams_transprobs(tprob_array, times = time_start, patient_wt = 1)$patient_wt,
               rep(1, prod(dim(tprob_array)[1:4])))
  expect_error(tparams_transprobs(tprob_array),
               paste0("'times' cannot be NULL if the number of time intervals ",
                      "is greater than 1"))
  expect_error(tparams_transprobs(tprob_array, times = time_start, 
                                  grp_id = rep(1, 3)),
                paste0("The length of 'grp_id' must be equal to the 3rd dimension of the ",
                       "array (i.e., the number of patients)."),
               fixed = TRUE)
  
})

test_that("tparams_transprobs returns array of matrices" , {
  expect_true(inherits(params_tprob$value, "array"))
  expect_equal(length(dim(params_tprob$value)), 3)
  expect_equal(dim(params_tprob$value)[1], dim(params_tprob$value)[2])
})

test_that("tparams_transprobs with only 1 time interval", {
  expect_true(inherits(params_tprob_t1, "tparams_transprobs"))
  expect_equal(params_tprob_t1$n_times, 1)
  expect_equal(nrow(params_tprob_t1$time_intervals), 1)
  expect_true(all(params_tprob_t1$time_id == 1))
})


# as.data.table ----------------------------------------------------------------
tprob_dt <- as.data.table(params_tprob)
  
test_that("as.data.table.tparams_transprobs returns a data.table" , {
  expect_true(inherits(as.data.table(params_tprob), "data.table"))
})

# tparams_transprobs.data.table------------------------------------------
params_tprob2 <- tparams_transprobs(tprob_dt)

test_that(paste0("tparams_transprobs returns the same values with ",
                 ".array and .data.table "), {
  expect_equal(params_tprob, params_tprob2)
  expect_equal(params_tprob_t1, 
               tparams_transprobs(tprob_dt[time_id == 1]))                
})

# Simulate model (from tparams object) -----------------------------------------
transmod <- CohortDtstmTrans$new(params = params_tprob)
econmod <- CohortDtstm$new(trans_model = transmod)

econmod$sim_stateprobs(n_cycles = 3)

test_that("CohortDtstmTrans$new automatically set 'start_stateprobs ",{
  tmp <- CohortDtstmTrans$new(params = params_tprob,
                              start_stateprobs = c(1, 0, 0))
  expect_equal(transmod$start_stateprobs, tmp$start_stateprobs)   
})

test_that("CohortDtstmTrans$sim_stateprobs() has correct grp_id ",{
  expect_true(all(econmod$stateprobs_$grp_id == 1))
})
                 
test_that("CohortDtstmTrans$sim_stateprobs() is correct ",{
  test_sim_stateprobs(econmod)
})

# Simulate model (from nnet object) --------------------------------------------
