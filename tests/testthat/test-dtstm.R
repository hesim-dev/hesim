context("dtstm.R unit tests")
library("data.table")
library("nnet")
library("msm")
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
tprob_array[, , 2, 1, 1, ] <- rdirichlet_mat(n_samples, make_alpha(tprob, .6))
tprob_array[, , 2, 1, 2, ] <- rdirichlet_mat(n_samples, make_alpha(tprob, .4))
tprob_array[, , 2, 2, 1, ] <- rdirichlet_mat(n_samples, make_alpha(tprob, .5))
tprob_array[, , 2, 2, 2, ] <- rdirichlet_mat(n_samples, make_alpha(tprob, .5))

tprob_array <- aperm(tprob_array, perm = c(6:3, 1, 2))

# tparams_transprobs.array() (6D) ----------------------------------------------
params_tprob <- tparams_transprobs(tprob_array, times = time_start)
params_tprob_t1 <- tparams_transprobs(tprob_array[,,,1,,, drop = FALSE])

test_that("Extra arguments with tparams_transprobs.array()" , {
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

test_that("tparams_transprobs() returns array of matrices" , {
  expect_true(inherits(params_tprob$value, "array"))
  expect_equal(length(dim(params_tprob$value)), 3)
  expect_equal(dim(params_tprob$value)[1], dim(params_tprob$value)[2])
})

test_that("tparams_transprobs() with only 1 time interval", {
  expect_true(inherits(params_tprob_t1, "tparams_transprobs"))
  expect_equal(params_tprob_t1$n_times, 1)
  expect_equal(nrow(params_tprob_t1$time_intervals), 1)
  expect_true(all(params_tprob_t1$time_id == 1))
})

# as.data.table.tparams_transprobs() -------------------------------------------
tprob_dt <- as.data.table(params_tprob)
  
test_that("as.data.table.tparams_transprobs() returns a data.table" , {
  expect_true(inherits(as.data.table(params_tprob), "data.table"))
})

# tparams_transprobs.data.table() ----------------------------------------------
params_tprob2 <- tparams_transprobs(tprob_dt)

test_that(paste0("tparams_transprobs() returns the same values with ",
                 ".array and .data.table "), {
  expect_equal(params_tprob, params_tprob2)
  expect_equal(params_tprob_t1, 
               tparams_transprobs(tprob_dt[time_id == 1]))                
})

test_that("tparams_transprobs.data.table() checks ID attributes", {
  expect_error(
    tparams_transprobs(tprob_dt[1:5]),
    paste0("The length of the ID variables is not consistent with the ",
            "number of unique values of each ID variable.")
  )
})

test_that("tparams_transprobs.data.table() throws error if there are no 'prob_' columns", {
  expect_error(
    tparams_transprobs(tprob_dt[, .(sample, strategy_id)]),
    "No columns with names starting with 'prob_'."
  )
})

# tparams_transprobs.tpmatrix() ------------------------------------------------
p <- c(.7, .6)
tpmat <- tpmatrix(
  C, p,
  0, 1
)
input_dat <- expand(hesim_dat)

test_that("tparams_transprobs() returns error if 'tpmatrix_id' has wrong class", {
  expect_error(
    tparams_transprobs(tpmat, 2),
    "'tpmatrix_id' must be of class 'data.frame'." 
  )
})

test_that("tparams_transprobs() returns error if 'tpmatrix_id' has wrong number of rows", {
  expect_error(
    tparams_transprobs(tpmat, data.frame(2)),
    "'object' and 'tpmatrix_id' must have the same number of rows." 
  )
})

test_that("tparams_transprobs() returns the correct class", {
  tpmat_id <- tpmatrix_id(input_dat, n_samples = 1)
  tp <- tpmatrix(
    C, c(.6, .7, .5, .4),
    0, 1
  )
  expect_true(
    inherits(tparams_transprobs(tp, tpmat_id),
             "tparams_transprobs")
  )
})

# tparams_transprobs.array() (3D) ----------------------------------------------
p <- c(.7, .6, .55, .58)
tpmat <- tpmatrix(
  C, p,
  0, 1
)
tparray <- as_array3(tpmat)
tpmat_id <- tpmatrix_id(input_dat, n_samples = 1)

test_that("tparams_transprobs.array() works as expected with 3D array", {
  tprob1 <- tparams_transprobs(tparray, tpmat_id)
  tprob2 <- tparams_transprobs(tpmat, tpmat_id)
  expect_equal(tprob1, tprob2)
})

test_that("tparams_transprobs.array() returns an error with the incorrect number of slices", {
  expect_error(tparams_transprobs(tparray[, , -1], tpmat_id),
               paste0("The third dimension of the array 'object' must equal ",
                      "the number or rows in 'tpmatrix_id'"))
})

# Initialize CohortDtstmTrans object -------------------------------------------
transmod <- CohortDtstmTrans$new(params = params_tprob)

test_that("CohortDtstmTrans$new() automatically sets 'start_stateprobs' ",{
  expect_equal(transmod$start_stateprobs, c(1, 0, 0))   
})

test_that("CohortDtstmTrans$new() 'start_stateprobs' normalizes to 1 ",{
  # Positive values
  v <- c(5, 5, 10, 10)
  tmp <- CohortDtstmTrans$new(params = params_tprob,
                              start_stateprobs = v)
  expect_equal(tmp$start_stateprobs, v/sum(v))
  tmp$start_stateprobs <- c(0, 0)
  expect_equal(tmp$start_stateprobs, c(1/2, 1/2))
  
  # All zeros
  tmp <- CohortDtstmTrans$new(params = params_tprob,
                              start_stateprobs = c(0, 0))
  expect_equal(tmp$start_stateprobs, c(1/2, 1/2))
})

test_that("CohortDtstmTrans$new() 'start_stateprobs' exceptions ",{
  expect_error(CohortDtstmTrans$new(params = params_tprob, 
                                    start_stateprobs = c(Inf, 1)),
               "Elements of 'state_stateprobs' cannot be infinite.")
  expect_error(CohortDtstmTrans$new(params = params_tprob, 
                                    start_stateprobs = c(0, -1)),
               "All elements of 'state_stateprobs' must be non-negative.")
})

test_that("CohortDtstmTrans 'trans_mat' must be a matrix of the correct form ",{
  msg_matrix <- "'trans_mat' must be a matrix"
  msg_form <- paste0("'trans_mat' is not of the correct form. Each row should ",
                     "contain integers from 0 to K - 1 where K is the number ", 
                     "of possible transitions (i.e., non-NA elements)")
  
  # Exceptions with $new()
  expect_error(CohortDtstmTrans$new(params = params_tprob, trans_mat = 1),
               msg_matrix)
  tmat_bad <- rbind(c(0, 0),
                    c(0, 0))
  expect_error(CohortDtstmTrans$new(params = params_tprob, trans_mat = tmat_bad),
              msg_form, fixed = TRUE)
  
  # Correct
  tmat_good <- rbind(c(0, 1, 2),
                     c(NA, 0, 1),
                     c(NA, NA, NA))
  tmp <- CohortDtstmTrans$new(params = params_tprob, trans_mat = tmat_good)
  expect_equal(tmp$trans_mat, tmat_good)
  
  # Active binding errors
  expect_error(tmp$trans_mat <- 1, msg_matrix)
  expect_error(tmp$trans_mat <- tmat_bad, msg_form, fixed = TRUE)
})

# Simulate model (from tparams object) -----------------------------------------
econmod <- CohortDtstm$new(trans_model = transmod)
econmod$sim_stateprobs(n_cycles = 3)

test_that("CohortDtstmTrans$sim_stateprobs() has correct grp_id ",{
  expect_true(all(econmod$stateprobs_$grp_id == 1))
})
                 
test_that("CohortDtstmTrans$sim_stateprobs() is correct ",{
  test_sim_stateprobs(econmod)
})

# Simulate model (from nnet object) --------------------------------------------
# Fit
transitions_data <- data.table(multinom3_exdata$transitions)
data_healthy <- transitions_data[state_from == "Healthy"]
fit_healthy <- multinom(state_to ~ strategy_name + female + age + year_cat,
                        data = data_healthy, trace = FALSE)
data_sick <- droplevels(transitions_data[state_from == "Sick"])
fit_sick <- multinom(state_to ~ strategy_name + female + age + year_cat,
                     data = data_sick, trace = FALSE)

# Construct model
## Setup
n_patients <- 100
patients <- transitions_data[year == 1, .(patient_id, age, female)][
  sample.int(nrow(transitions_data[year == 1]), n_patients)][
    , grp_id := 1:n_patients]
hesim_dat <- hesim_data(
  patients = patients,
  strategies = data.table(strategy_id = 1:2,
                          strategy_name = c("Reference", "Intervention")),
  states = data.table(state_id = c(1, 2),
                      state_name = c("Healthy", "Sick")) # Non-death health states
)

## The model
n_samples <- 10
tmat <- rbind(c(0, 1, 2),
              c(NA, 0, 1),
              c(NA, NA, NA))
transfits <- multinom_list(healthy = fit_healthy, sick = fit_sick)
tintervals <- time_intervals(unique(transitions_data[, .(year_cat)])
                             [, time_start := c(0, 2, 6)])
transmod_data <- expand(hesim_dat, times = tintervals)
transmod <- create_CohortDtstmTrans(transfits,
                                    input_data = transmod_data,
                                    trans_mat = tmat,
                                    n = n_samples,
                                    uncertainty = "none")


test_that(paste0("create_CohortDtstmTrans$sim_stateprobs() is consistent with ",
                 "predict.multinom()"), {
  hesim_probs <- transmod$sim_stateprobs(n_cycles = 1)[t == 1]
  hesim_probs[, state_id := factor(state_id, labels = c("Healthy", "Sick", "Dead"))]
  hesim_probs <- dcast(hesim_probs, 
                       strategy_id + patient_id ~ state_id,
                       value.var = "prob")          
  multinom_probs <- predict(fit_healthy, newdata = transmod_data[time_id == 1],
                           type = "prob")
  rownames(multinom_probs) <- NULL
  expect_equal(multinom_probs,
               as.matrix(hesim_probs[, c("Healthy", "Sick", "Dead")]))
})
     
test_that(paste0("create_CohortDtstmTrans$sim_stateprobs() with mulinom() objects"), {
  tpmatrix_multinom <- function(fits, data, patid){
    newdata <- data[patient_id == patid]
    n_times <- length(unique(transmod_data$time_id))
    tpmatrix1 <- tpmatrix2 <- array(NA, dim = c(3, 3, n_times))
    for (j in 1:n_times){
      probs_healthy <- predict(fits$healthy, newdata[time_id == j], type = "probs")
      probs_sick <- predict(fits$sick, newdata[time_id == j], type = "probs")
      tpmatrix1[, , j] <- rbind(probs_healthy[1, ],
                                c(0, 1 - probs_sick[1], probs_sick[1]),
                                c(0, 0, 1))
      tpmatrix2[, , j] <- rbind(probs_healthy[2, ],
                                c(0, 1 - probs_sick[2], probs_sick[2]),
                                c(0, 0, 1))
    }
    return(list(p_ref = tpmatrix1,
                p_int = tpmatrix2))  
  } 
  times <- 0:10
  patid <- sample(unique(transmod_data$patient_id), 1)
  p <- tpmatrix_multinom(transfits, transmod_data, patid = patid)
  stprobs_ref <- sim_markov_chain(p = p$p_ref,
                                  x0 = c(1, 0, 0),
                                  times = times,
                                  time_stop = unique(transmod_data$time_stop))
  stprobs_int <- sim_markov_chain(p = p$p_int,
                                  x0 = c(1, 0, 0),
                                  times = times,
                                  time_stop = unique(transmod_data$time_stop))
  
  hesim_stprobs <- transmod$sim_stateprobs(n_cycles = max(times))[patient_id == patid]
  hesim_stprobs <- dcast(hesim_stprobs, 
                         strategy_id + patient_id + t~ state_id,
                         value.var = "prob") 
  test_equal <- function(R_stprobs, hesim_stprobs, strat_id){
    hesim_stprobs <- as.matrix(hesim_stprobs[strategy_id == strat_id][,
                                  c("1", "2", "3"), with = FALSE])
    colnames(hesim_stprobs) <- NULL
    rownames(hesim_stprobs) <- times
    expect_equal(hesim_stprobs, R_stprobs)
  }
  test_equal(stprobs_ref, hesim_stprobs, strat_id = 1)
  test_equal(stprobs_int, hesim_stprobs, strat_id = 2)
})

test_that(paste0("create_CohortDtstmTrans does not support offset term ", 
                 "with mulinom() objects"), {
  m <- matrix(c(1, 100, 1), nrow = nrow(data_healthy), ncol = 3, byrow = TRUE)
  fit_healthy2 <- multinom(state_to ~ offset(m) + strategy_name + female + age +
                             year_cat,
                          data = data_healthy, trace = FALSE)  
  transfits2 <- multinom_list(healthy = fit_healthy2, sick = fit_sick)
  expect_error(create_CohortDtstmTrans(transfits2,
                                       input_data = transmod_data,
                                       trans_mat = tmat),
               "An offset is not supported")
})

# Create model from a msm object -----------------------------------------------
set.seed(101)
strategies <- data.table(strategy_id = c(1, 2, 3),
                         strategy_name = factor(c("SOC", "New 1", "New 2")))
patients <- data.table(patient_id = 1:2)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)
transmod_data <- expand(hesim_dat)
qinit <- rbind(
  c(0, 0.28163, 0.01239),
  c(0, 0, 0.10204),
  c(0, 0, 0)
)
fit <- msm(state_id ~ time, subject = patient_id, 
           data = onc3p[patient_id %in% sample(patient_id, 100)],
           covariates = list("1-2" =~ strategy_name), 
           qmatrix = qinit)

test_that("create_CohortDtstmTrans.msm() returns correct transition probability matrices with no uncertainty", {
 
  transmod <- create_CohortDtstmTrans(fit,
                                      input_data = transmod_data,
                                      cycle_length = 1/2,
                                      fixedpars = 2,
                                      uncertainty = "none")
  expect_equal(transmod$params$n_samples, 1)
  expect_equal(
    expmat(qmatrix(fit, transmod_data, uncertainty = "none"), t = 1/2),
    transmod$params$value
  )
})

test_that(paste0("create_CohortDtstmTrans.msm() returns transition probability matrices",
                 "with correct dimensions when there is uncertainty"), {

  transmod <- create_CohortDtstmTrans(fit,
                                      input_data = transmod_data,
                                      cycle_length = 1/2,
                                      fixedpars = 2,
                                      n = 2)
  expect_equal(transmod$params$n_samples, 2)
  expect_equal(dim(transmod$params$value)[3], 2 * nrow(transmod_data))
})


# Create model from a model_def object -----------------------------------------
test_that(paste0("define_model() works to create CohortDtstmTrans object"), {
  hesim_dat <- hesim_data(
    strategies = data.table(strategy_id = 1:2,
                            strategy_name = c("Monotherapy", "Combination therapy")),
    patients = data.table(patient_id = 1)
  )
  data <- expand(hesim_dat)
  
  # Define the model
  rng_def <- define_rng({
    alpha <- matrix(c(1251, 350, 116, 17,
                      0, 731, 512, 15,
                      0, 0, 1312, 437,
                      0, 0, 0, 469),
                    nrow = 4, byrow = TRUE)
    rownames(alpha) <- colnames(alpha) <- c("A", "B", "C", "D")
    lrr_mean <- log(.509)
    lrr_se <- (log(.710) - log(.365))/(2 * qnorm(.975))
    list(
      p_mono = dirichlet_rng(alpha),
      rr_comb = lognormal_rng(lrr_mean, lrr_se),
      u = 1,
      c_zido = 2278,
      c_lam = 2086.50,
      c_med = gamma_rng(mean = c(A = 2756, B = 3052, C = 9007),
                        sd = c(A = 2756, B = 3052, C = 9007))
    )
  }, n = 2)
  
  tparams_def <- define_tparams({
    rr = ifelse(strategy_name == "Monotherapy", 1, rr_comb)
    list(
      tpmatrix = tpmatrix(
        C, p_mono$A_B * rr, p_mono$A_C * rr, p_mono$A_D * rr,
        0, C, p_mono$B_C * rr, p_mono$B_D * rr,
        0, 0, C, p_mono$C_D * rr,
        0, 0, 0, 1),
      utility = u,
      costs = list(
        drug = ifelse(strategy_name == "Monotherapy",
                      c_zido, c_zido + c_lam),
        medical = c_med
      ) 
    )
  })
  
  model_def <- define_model(
    tparams_def = tparams_def,
    rng_def = rng_def)
  econmod <- create_CohortDtstm(model_def, data)
  
  # Test
  expect_true(inherits(econmod, "CohortDtstm"))
})

