context("input_data unit tests")
library("flexsurv")
library("data.table")
rm(list = ls())

dt_strategies <- data.table(strategy_id = c(1, 2))
dt_patients <- data.table(patient_id = seq(1, 3), 
                          age = c(45, 47, 60),
                          female = c(1, 0, 0),
                          group = factor(c("Good", "Medium", "Poor")))
dt_lines <- create_lines_dt(list(c(1, 2, 5), c(1, 2)))
dt_states <- data.frame(state_id =  seq(1, 3),
                        state_name = factor(paste0("state", seq(1, 3))))
dt_trans <- data.frame(transition_id = seq(1, 4),
                       from = c(1, 1, 2, 2),
                       to = c(2, 3, 1, 3))
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients,
                        lines = dt_lines,
                        states = dt_states,
                        transitions = dt_trans)

# create_lines_dt --------------------------------------------------------------
test_that("create_lines_dt", {
  dt_lines <- create_lines_dt(list(c(1, 2, 5), c(1, 2)))
  
  expect_true(inherits(dt_lines, "data.table"))
  expect_equal(dt_lines$treatment_id[3], 5)
  expect_equal(dt_lines$line, 
               c(seq(1, 3), seq(1, 2)))
  
  # explicit strategy ids
  dt_lines <- create_lines_dt(list(c(1, 2, 5), c(1, 2)),
                              strategy_ids = c(3, 5))
  expect_equal(dt_lines$strategy_id, c(3, 3, 3, 5, 5))
  
  # errors
  expect_error(input_data$lines_dt(list(c("tx1", "tx2"),
                                  c("tx1"))))
})

# create_trans_dt --------------------------------------------------------------
test_that("create_trans_dt", {
  tmat <- rbind(c(NA, 1, 2),
                c(NA, NA, 3),
                c(NA, NA, NA))
  dt_trans <- create_trans_dt(tmat)
  
  expect_true(inherits(dt_trans, "data.table"))
  expect_equal(dt_trans$transition_id, 
               c(1, 2, 3))
  expect_equal(dt_trans$from, 
               c(1, 1, 2))
  expect_equal(dt_trans$to, 
               c(2, 3, 3))
  
  # Row and column names
  rownames(tmat) <- c("No BOS", "BOS", "Death")
  dt_trans <- create_trans_dt(tmat)
  expect_equal(dt_trans$from_name, NULL)
  
  colnames(tmat) <- rownames(tmat)
  dt_trans <- create_trans_dt(tmat)
  expect_equal(dt_trans$from_name, rownames(tmat)[c(1, 1, 2)])
  expect_equal(dt_trans$to_name, colnames(tmat)[c(2, 3, 3)])
})

# hesim data -------------------------------------------------------------------
test_that("hesim_data", {

  # strategy by patient
  hesim_dat <- hesim_data(strategies = dt_strategies,
                          patients = dt_patients)
  
  expect_true(inherits(hesim_dat, "hesim_data"))
  expect_equal(hesim_dat$state, NULL)
  expect_equal(hesim_dat$patients, dt_patients)
  
  # strategy by patient by state
  hesim_dat <- hesim_data(strategies = dt_strategies,
                          patients = dt_patients, 
                          states = dt_states)
  expect_equal(hesim_dat$states, dt_states)
  
  # Expand
  expanded_dt <- expand_hesim_data(hesim_dat, by = c("strategies"))
  expect_equal(expanded_dt$data, data.table(dt_strategies))
  expect_equal(expanded_dt$id_vars, "strategy_id")
  expanded_dt <- expand_hesim_data(hesim_dat, by = c("strategies", "patients"))
  expanded_dt2 <- expand_hesim_data(hesim_dat, by = c("patients", "strategies"))
  expect_equal(nrow(expanded_dt$data), 
               nrow(dt_strategies) * nrow(dt_patients))
  expect_equal(expanded_dt$data, expanded_dt2$data)
  expect_equal(expanded_dt$id_vars, expanded_dt2$id_vars)
  expect_equal(expanded_dt$id_vars, c("strategy_id", "patient_id"))
  
  # errors
  expect_error(expand_hesim_data(hesim_dat, by = c("strategies", "patients", 
                                                  "states", "transitions")))
  expect_error(expand_hesim_data(hesim_dat, by = c("strategies", "patients", 
                                                  "states", "wrong_table")))
  hesim_dat2 <- hesim_dat[c("strategies", "patients")]
  expect_error(expand_hesim_data(hesim_dat2, by = c("strategies", "patients", 
                                                  "states")))
})

# input_data class -------------------------------------------------------------
# By treatment strategy and patient
dat <- expand_hesim_data(hesim_dat)$data
input_dat <- input_data(X = list(mu = model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)))

## X must be a list
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id))))

## X must be a list of matrices
expect_error(input_data(X = list(model.matrix(~ age, dat), 2),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id))))

## Number of rows in X is inconsistent with strategy_id 
expect_error(input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id[-1],
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id))))

## Size of patient_id is incorrect
expect_error(input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = sort(dat$patient_id),
                       n_patients = length(unique(dat$strategy_id))))

## n_patients is incorrect v1
expect_error(input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$strategy_id,
                       n_patients = length(unique(dat$strategy_id))))

## n_patients is incorrect v2
expect_error(input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = 1))

## patient_id is not sorted correctly
expect_error(input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = sort(dat$patient_id),
                       n_patients = length(unique(dat$patient_id))))


# By treatment strategy, line, and patient
dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "lines"))$data
n_lines <- hesim_dat$lines[, .N, by = "strategy_id"]
input_dat <- input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)),
                       line = dat$line,
                       n_lines = n_lines)

## n_lines is incorrect v1
n_lines[, N := N + 1]
expect_error(input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)),
                       line = dat$line,
                       n_lines = n_lines))

## n_lines is incorrect v2
n_lines[, N := N - 1]
expect_error(input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)),
                       line = dat$line,
                       n_lines = n_lines + 1))

## line is not sorted correctly
expect_error(input_data(X = list(model.matrix(~ age, dat)),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)),
                       line = sort(dat$line),
                       n_lines = n_lines))

# create_input_data with formula objects ---------------------------------------
test_that("create_input_data.formula_list", {
  dat <- expand_hesim_data(hesim_dat)
  f_list <- formula_list(list(f1 = formula(~ age), f2 = formula(~ 1)))
  expect_equal(class(f_list), "formula_list")
  input_dat <- create_input_data(f_list, dat)
  
  expect_equal(length(input_dat$X), length(f_list))
  expect_equal(names(input_dat$X), names(f_list))
  expect_equal(as.numeric(input_dat$X$f1[, "age"]), dat$data$age)
  expect_equal(ncol(input_dat$X$f1), 2)
  expect_equal(ncol(input_dat$X$f2), 1)
})

# create_input_data with lm objects --------------------------------------------
dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "states"))
fit1 <- stats::lm(costs ~ female + state_name, data = psm4_exdata$costs$medical)

test_that("create_input_data.lm", {
  input_dat <- create_input_data(fit1, dat)
  
  expect_equal(ncol(input_dat$X$mu), 4)
  expect_equal(as.numeric(input_dat$X$mu[, "female"]), dat$data$female)
})

test_that("create_input_data.lm_list", {
  fit2 <- stats::lm(costs ~ 1, data = psm4_exdata$costs$medical)
  fit_list <- hesim:::lm_list(fit1 = fit1, fit2 = fit2)
  input_dat <- create_input_data(fit_list, dat)
  
  expect_equal(ncol(input_dat$X$fit1$mu), 4)
  expect_equal(ncol(input_dat$X$fit2$mu), 1)
  expect_equal(as.numeric(input_dat$X$fit1$mu[, "female"]), dat$data$female)
})

# create_input_data with flexsurvreg objects -----------------------------------
test_that("create_input_data.flexsurv", {
  dat <- expand_hesim_data(hesim_dat)
  fit <- flexsurv::flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                              anc = list(sigma = ~ group), 
                              dist = "gengamma") 
  input_dat <- create_input_data(fit, dat)
  
  expect_equal(input_dat$strategy_id, dat$data$strategy_id)
  expect_equal(input_dat$state_id, dat$data$state_id)
  expect_equal(input_dat$patient_id, dat$data$patient_id)
  expect_equal(class(input_dat$X), "list")
  expect_equal(class(input_dat$X[[1]]), "matrix")
  expect_equal(length(input_dat$X), 3)
  expect_equal(ncol(input_dat$X$mu), 3)
  expect_equal(ncol(input_dat$X$sigma), 3)
  expect_equal(ncol(input_dat$X$Q), 1)
})

fit1_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                                  data = ovarian, dist = "weibull")
fit1_exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                                  data = ovarian, dist = "exp")
flexsurvreg_list1 <- flexsurvreg_list(wei = fit1_wei, exp = fit1_exp)
dat <- expand_hesim_data(hesim_dat)

test_that("create_input_data.flexsurv_list", {
  input_dat <- create_input_data(flexsurvreg_list1, dat)  
  
  expect_equal(class(input_dat$X$wei$shape), "matrix")
})

fit2_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1 + age, 
                                  data = ovarian, 
                                  dist = "weibull")
fit2_exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1 + age, 
                                  data = ovarian, 
                                  dist = "exp")
flexsurvreg_list2 <- flexsurvreg_list(wei = fit2_wei, exp = fit2_exp)
joined_flexsurvreg_list <- joined_flexsurvreg_list(mod1 = flexsurvreg_list1,
                                                   mod2 = flexsurvreg_list2,
                                                   times = list(2, 5))

test_that("create_input_data.joined_flexsurv_list", {
  input_dat <- create_input_data(joined_flexsurvreg_list, dat)  
  
  expect_equal(input_dat$state_id, dat$state_id)
  expect_equal(class(input_dat$X[[1]]$wei$shape), "matrix")
})
