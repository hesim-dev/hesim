context("input_data unit tests")
library("flexsurv")
library("data.table")
rm(list = ls())

strategies <- data.table(strategy_id = c(1, 2))
patients <- data.table(
  patient_id = seq(1, 3), 
  age = c(45, 47, 60),
  female = c(1, 0, 0),
  group = factor(c("Good", "Medium", "Poor")),
  ecog = c("Asymptomatic", "Symptomatic (ambulatory)",
           "In bed <50%")
)
states <- data.table(
  state_id =  seq(1, 3),
  state_name = factor(paste0("state", seq(1, 3)))
)
transitions <- data.table(
  transition_id = seq(1, 4),
  from = c(1, 1, 2, 2),
  to = c(2, 3, 1, 3)
)
hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients,
  states = states,
  transitions = transitions
)
input_data <- expand(hesim_dat)

# input_mats class works as expected -------------------------------------------
im <- input_mats(
  X = list(mu = model.matrix(~ age, input_data)),
  strategy_id = input_data$strategy_id,
  n_strategies = length(unique(input_data$strategy_id)),
  patient_id = input_data$patient_id,
  n_patients = length(unique(input_data$patient_id))
)

test_that("input_mats() works as expected", {
  expected_X <- as.matrix(data.frame(1, input_data$age))
  colnames(expected_X) <- c("(Intercept)", "age")
  
  expect_true(inherits(im, "input_mats"))
  expect_equivalent(im$X$mu, expected_X)
})

test_that("print.input_mats() works as expected", {
  expect_output(print(im), "An \"input_mats\" object")
  expect_output(print(im), paste0("Column binding the ID variables with all ",
                                  "variables contained in the X matrices:"))
  expect_output(print(im), "Number of unique values of ID variables:")
})

test_that("as.data.table.input_mats() works as expected", {
 imd <- as.data.table(im)
 expect_true(inherits(imd, "data.table"))
 expect_equal(imd$age, input_data$age)
})

# input_mats class throws errors -----------------------------------------------
test_that("input_mats() throws error if X is not a list", {
  expect_error(
    input_mats(
      X = model.matrix(~ age, input_data),
      strategy_id = input_data$strategy_id,
      n_strategies = length(unique(input_data$strategy_id)),
      patient_id = input_data$patient_id,
      n_patients = length(unique(input_data$patient_id))
    ),
    "'X' must be a list or a list of lists."
  )
})

test_that("input_mats() throws error if X is not a list of matrices", {
  expect_error(
    input_mats(
      X = list(2),
      strategy_id = input_data$strategy_id,
      n_strategies = length(unique(input_data$strategy_id)),
      patient_id = input_data$patient_id,
      n_patients = length(unique(input_data$patient_id))
    ),
    "'X' must be a list or list of lists of matrices."
  )
})

test_that("input_mats() throws error if the rows in X are inconsistent with strategy_id", {
  expect_error(
    input_mats(
      X = list(mu = model.matrix(~ age, input_data)),
      strategy_id = input_data$strategy_id[-1],
      n_strategies = length(unique(input_data$strategy_id)),
      patient_id = input_data$patient_id,
      n_patients = length(unique(input_data$patient_id))
    ),
    "The length of 'strategy_id' does not equal the number of rows in the 'X' matrices."
  )
})

# create_input_mats with lm objects --------------------------------------------
input_data <- expand(hesim_dat, by = c("strategies", "patients", "states"))
fit1 <- lm(costs ~ female + state_name, data = psm4_exdata$costs$medical)

test_that("create_input_mats.lm() works with both data.table and data.frame input data", {
  # With data.table input data
  im1 <- create_input_mats(fit1, input_data)
  expect_equal(ncol(im1$X$mu), 4)
  expect_equal(as.numeric(im1$X$mu[, "female"]), input_data$female)
  
  # With data.frame input data
  input_data2 <- copy(input_data)
  setattr(input_data2, "class", c("expanded_hesim_data", "data.frame"))
  im2 <- create_input_mats(fit1, input_data2)
  expect_equal(im1, im2)
})

test_that("create_input_mats.lm() works with times", {
  d <- expand(hesim_dat, by = c("strategies", "patients", "states"), 
              times = c(0, 2))
  im <- create_input_mats(fit1, d)
  
  expect_output(print(im), "Time intervals:")
  imd <- as.data.table(im)
  expect_true(
    all(c("strategy_id", "patient_id", "state_id", "time_id", "time_start", 
          "time_stop", "state_namestate2", "state_namestate3") %in%
        colnames(imd))
  )
  expect_equal(d$female, imd$female)
})

test_that("create_input_mats.lm() thows error if input data is not a data.table or data.frame", {
  d <- copy(input_data)
  setattr(d, "class", "expanded_hesim_data")
  expect_error(
    create_input_mats(fit1, d),
    "'input_data' must inherit from either 'data.table' or 'data.frame'"
  )
})

# create_input_mats with params_lm objects -------------------------------------
test_that("create_input_mats.params_lm() works as expected", {
  p <- params_lm(
    coef = data.frame(intercept = c(.2, .3), age = c(.02, .05))
  )
  id <- expand(hesim_dat)[, intercept := 1]
  im <- create_input_mats(p, id)
  expect_equal(im$X$mu[, "age"], id$age)
  expect_equal(im$patient_id, id$patient_id)
})

# create_input_mats with flexsurvreg objects -----------------------------------
input_data <- expand(hesim_dat)[, intercept := 1]

test_that("input_mats.flexsurvreg() returns the correct columns", {
  lung <- data.table(survival::lung)
  lung[, status := ifelse(status == 2, 0, 1)]
  lung[, ecog := factor(
    ph.ecog, levels = 0:3,
    labels = c("Asymptomatic", "Symptomatic (ambulatory)",
               "In bed <50%", "In bed > 50% of the day")
  )]
  fit <- flexsurvreg(Surv(time, status) ~ poly(age, 2) + factor(ecog), 
                     data = lung,
                     dist = "weibull") 
  p <- create_params(fit)
  terms <- colnames(p$coefs$scale)
  im <- create_input_mats(fit, input_data)
  expect_equal(colnames(im$X$scale)[-1], terms[-1])
})

test_that(paste0("create_input_mats.flexsurvreg() works with regression ",
                  "coefficients on ancillary parameters"), {
  fit <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                    anc = list(sigma = ~ group), 
                    dist = "gengamma") 
  im <- create_input_mats(fit, input_data)
  
  expect_equal(im$strategy_id, input_data$strategy_id)
  expect_equal(im$state_id, input_data$state_id)
  expect_equal(im$patient_id, input_data$patient_id)
  expect_equal(class(im$X), "list")
  expect_true(inherits(im$X[[1]], "matrix"))
  expect_equal(length(im$X), 3)
  expect_equal(ncol(im$X$mu), 3)
  expect_equal(ncol(im$X$sigma), 3)
  expect_equal(ncol(im$X$Q), 1)
})

test_that("create_input_mats.flexsurv_list() works as expected", {
  fw <- flexsurvreg(Surv(futime, fustat) ~ 1, 
                    data = ovarian, dist = "weibull")
  fe <- flexsurvreg(Surv(futime, fustat) ~ 1, 
                    data = ovarian, dist = "exp")
  fl <- flexsurvreg_list(wei = fw, exp = fe)
  
  im <- create_input_mats(fl, input_data)  
  expect_true(inherits(im, "input_mats"))
  expect_true(inherits(im$X$wei$shape, "matrix"))
})

# create_input_mats with params_surv objects -----------------------------------
p_wei <- params_surv(
  coef = list(
    scale = data.frame(intercept = c(.2, .3), 
                       age = c(.02, .05)),
    shape = data.frame(intercept = c(.2, .3))
  ),
  dist = "weibull"
) 

test_that("create_input_mats.params_surv() works as expected", {
  im <- create_input_mats(p_wei, input_data)
  expect_equal(im$X$shape[, "intercept"], input_data$intercept)
  expect_equal(im$X$scale[, "age"], input_data$age)
  expect_equal(im$strategy_id, input_data$strategy_id)
})

test_that("create_input_mats.params_surv_list() works as expected", {
  p_exp <- params_surv(
    coef = list(
      rate = data.frame(intercept = c(.2, .3), 
                         age = c(.02, .05))
    ),
    dist = "exp"
  ) 
  p <- params_surv_list(wei = p_wei, exp = p_exp)
  im <- create_input_mats(p, input_data) 
  expect_equal(im$X$wei$scale[, "age"], input_data$age)
  expect_equal(im$X$exp$rate[, "age"], input_data$age)
  expect_equal(im$strategy_id, input_data$strategy_id)
})


