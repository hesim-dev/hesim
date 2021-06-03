context("input_data unit tests")
library("flexsurv")
library("data.table")
rm(list = ls())

strategies_dt <- data.table(strategy_id = c(1, 2))
patients_dt <- data.table(
  patient_id = seq(1, 3), 
  age = c(45, 47, 60),
  female = c(1, 0, 0),
  group = factor(c("Good", "Medium", "Poor")),
  ecog = c("Asymptomatic", "Symptomatic (ambulatory)",
           "In bed <50%")
)
states_dt <- data.frame(
  state_id =  seq(1, 3),
  state_name = factor(paste0("state", seq(1, 3)))
)
trans_dt <- data.frame(
  transition_id = seq(1, 4),
  from = c(1, 1, 2, 2),
  to = c(2, 3, 1, 3)
)
hesim_dat <- hesim_data(
  strategies = strategies_dt,
  patients = patients_dt,
  states = states_dt,
  transitions = trans_dt
)

# input_mats class -------------------------------------------------------------
# By treatment strategy and patient
dat <- expand(hesim_dat)
X <- input_mats(X = list(mu = model.matrix(~ age, dat)),
                strategy_id = dat$strategy_id,
                n_strategies = length(unique(dat$strategy_id)),
                patient_id = dat$patient_id,
                n_patients = length(unique(dat$patient_id)))

## X must be a list
expect_error(input_mats(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id))))

## X must be a list of matrices
expect_error(input_mats(X = list(model.matrix(~ age, dat), 2),
                        strategy_id = dat$strategy_id,
                        n_strategies = length(unique(dat$strategy_id)),
                        patient_id = dat$patient_id,
                        n_patients = length(unique(dat$patient_id))))

## Number of rows in X is inconsistent with strategy_id 
expect_error(input_mats(X = list(model.matrix(~ age, dat)),
                        strategy_id = dat$strategy_id[-1],
                        n_strategies = length(unique(dat$strategy_id)),
                        patient_id = dat$patient_id,
                        n_patients = length(unique(dat$patient_id))))

# create_input_mats with lm objects or params_lm objects -----------------------
dat <- expand(hesim_dat, by = c("strategies", "patients", "states"))
fit1 <- stats::lm(costs ~ female + state_name, data = psm4_exdata$costs$medical)

test_that("create_input_mats.lm", {
  input_mats1 <- create_input_mats(fit1, dat)
  expect_equal(ncol(input_mats1$X$mu), 4)
  expect_equal(as.numeric(input_mats1$X$mu[, "female"]), dat$female)
  
  # Works with data.frame
  dat_df = copy(dat)
  setattr(dat_df, "class", c("expanded_hesim_data", "data.frame"))
  input_mats2 <- create_input_mats(fit1, dat_df)
  expect_equal(input_mats1, input_mats2)
  
  # Error if not data.table or data.frame
  setattr(dat_df, "class", "expanded_hesim_data")
  expect_error(create_input_mats(fit1, dat_df))
})

test_that("create_input_mats.lm_list", {
  fit2 <- stats::lm(costs ~ 1, data = psm4_exdata$costs$medical)
  fit_list <- hesim:::lm_list(fit1 = fit1, fit2 = fit2)
  input_mats <- create_input_mats(fit_list, dat)
  
  expect_equal(ncol(input_mats$X$fit1$mu), 4)
  expect_equal(ncol(input_mats$X$fit2$mu), 1)
  expect_equal(as.numeric(input_mats$X$fit1$mu[, "female"]), dat$female)
})

test_that("create_input_mats.params_lm", {
  coef <- as.matrix(data.frame(intercept = c(.2, .3), age = c(.02, .05)))
  params <- params_lm(coef = coef)
  data <- data.table(intercept = c(1, 1), age = c(55, 65),
                     patient_id = c(1, 2), strategy_id = c(1, 1))
  setattr(data, "id_vars", c("patient_id", "strategy_id"))
  setattr(data, "class", c("expanded_hesim_data", "data.table", "data.frame"))
  input_mats <- create_input_mats(params, data)
  expect_equal(input_mats$X$mu[, "intercept"], c(1, 1))
  expect_equal(input_mats$patient_id, c(1, 2))
})

# create_input_mats with flexsurvreg or params_surv objects --------------------
dat <- expand(hesim_dat)

test_that("input_mats.flexsurvreg returns the correct columns", {
  lung <- data.table(survival::lung)
  lung[, status := ifelse(status == 2, 0, 1)]
  lung[, ecog := factor(
    ph.ecog, levels = 0:3,
    labels = c("Asymptomatic", "Symptomatic (ambulatory)",
               "In bed <50%", "In bed > 50% of the day")
  )]
  fit <- flexsurv::flexsurvreg(Surv(time, status) ~ poly(age, 2) + factor(ecog), 
                               data = lung,
                               dist = "weibull") 
  params <- create_params(fit)
  terms <- colnames(params$coefs$scale)
  input_mats <- create_input_mats(fit, dat)
  expect_equal(colnames(input_mats$X$scale)[-1], terms[-1])
})

test_that(paste0("create_input_mats.flexsurvreg works with regression ",
                  "coefficients on ancillary parameters"), {
  fit <- flexsurv::flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                              anc = list(sigma = ~ group), 
                              dist = "gengamma") 
  input_mats <- create_input_mats(fit, dat)
  
  expect_equal(input_mats$strategy_id, dat$strategy_id)
  expect_equal(input_mats$state_id, dat$state_id)
  expect_equal(input_mats$patient_id, dat$patient_id)
  expect_equal(class(input_mats$X), "list")
  expect_true(inherits(input_mats$X[[1]], "matrix"))
  expect_equal(length(input_mats$X), 3)
  expect_equal(ncol(input_mats$X$mu), 3)
  expect_equal(ncol(input_mats$X$sigma), 3)
  expect_equal(ncol(input_mats$X$Q), 1)
})

fit1_wei <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ 1, 
                                  data = ovarian, dist = "weibull")
fit1_exp <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ 1, 
                                  data = ovarian, dist = "exp")
flexsurvreg_list1 <- flexsurvreg_list(wei = fit1_wei, exp = fit1_exp)

test_that("create_input_mats.flexsurv_list", {
  input_mats <- create_input_mats(flexsurvreg_list1, dat)  
  expect_true(inherits(input_mats$X$wei$shape, "matrix"))
})

test_that("create_input_mats.params_surv", {
  # params_surv
  coef_wei <- list(scale = as.matrix(data.frame(intercept = c(.2, .3), 
                                            age = c(.02, .05))),
               shape = as.matrix(data.frame(intercept = c(.2, .3))))
  params_wei <- params_surv(coef = coef_wei,
                        dist = "weibull") 
  data <- data.table(intercept = c(1, 1), age = c(55, 65),
                     patient_id = c(1, 2), strategy_id = c(1, 1))
  setattr(data, "id_vars", c("patient_id", "strategy_id"))
  setattr(data, "class", c("expanded_hesim_data", "data.table", "data.frame"))  
  input_mats <- create_input_mats(params_wei, data)
  expect_equal(input_mats$X$shape[, "intercept"], c(1, 1))
  expect_equal(input_mats$X$scale[, "age"], data$age)
  expect_equal(input_mats$strategy_id, data$strategy_id)
  
  # params_surv_list
  coef_exp <- list(rate = as.matrix(data.frame(intercept = c(.2, .3), 
                                            age = c(.02, .05))))
  params_exp <- params_surv(coef = coef_exp,
                                dist = "exp") 
  params <- params_surv_list(wei = params_wei, exp = params_exp)
  input_mats <- create_input_mats(params, data) 
  expect_equal(input_mats$X$wei$scale[, "age"], data$age)
  expect_equal(input_mats$X$exp$rate[, "age"], data$age)
  expect_equal(input_mats$strategy_id, data$strategy_id)
})

