context("input_data unit tests")
library("flexsurv")
library("data.table")

dt.strategies <- data.table(strategy_id = c(1, 2))
dt.patients <- data.table(patient_id = seq(1, 3), 
                          age = c(45, 47, 60),
                          female = c(1, 0, 0),
                          group = factor(c("Good", "Medium", "Poor")))
dt.lines <- lines_dt(list(c(1, 2, 5), c(1, 2)))
dt.states <- data.frame(state_id =  seq(1, 3),
                         state_name = factor(paste0("state", seq(1, 3))))
dt.trans <- data.frame(transition_id = seq(1, 4),
                       from = c(1, 1, 2, 2),
                       to = c(2, 3, 1, 3))
hesim.dat <- hesim_data(strategies = dt.strategies,
                        patients = dt.patients,
                        lines = dt.lines,
                        states = dt.states,
                        transitions = dt.trans)

# lines_dt ---------------------------------------------------------------------
test_that("lines_dt", {
  dt.lines <- lines_dt(list(c(1, 2, 5), c(1, 2)))
  
  expect_true(inherits(dt.lines, "data.table"))
  expect_equal(dt.lines$treatment_id[3], 5)
  expect_equal(dt.lines$line, 
               c(seq(1, 3), seq(1, 2)))
  
  # explicit strategy ids
  dt.lines <- lines_dt(list(c(1, 2, 5), c(1, 2)),
                                  strategy_ids = c(3, 5))
  expect_equal(dt.lines$strategy_id, c(3, 3, 3, 5, 5))
  
  # errors
  expect_error(input_data$lines_dt(list(c("tx1", "tx2"),
                                  c("tx1"))))
})

# hesim data -------------------------------------------------------------------
test_that("hesim_data", {

  # strategy by patient
  hesim.dat <- hesim_data(strategies = dt.strategies,
                       patients = dt.patients)
  
  expect_true(inherits(hesim.dat, "hesim_data"))
  expect_equal(hesim.dat$state, NULL)
  expect_equal(hesim.dat$patients, dt.patients)
  
  # strategy by patient by state
  hesim.dat <- hesim_data(strategies = dt.strategies,
                       patients = dt.patients, 
                       states = dt.states)
  expect_equal(hesim.dat$states, dt.states)
  
  # Expand
  expanded.dt <- expand_hesim_data(hesim.dat, by = c("strategies"))
  expect_equal(expanded.dt$data, data.table(dt.strategies))
  expect_equal(expanded.dt$id_vars, "strategy_id")
  expanded.dt <- expand_hesim_data(hesim.dat, by = c("strategies", "patients"))
  expanded.dt2 <- expand_hesim_data(hesim.dat, by = c("patients", "strategies"))
  expect_equal(nrow(expanded.dt$data), 
               nrow(dt.strategies) * nrow(dt.patients))
  expect_equal(expanded.dt$data, expanded.dt2$data)
  expect_equal(expanded.dt$id_vars, expanded.dt2$id_vars)
  expect_equal(expanded.dt$id_vars, c("strategy_id", "patient_id"))
  
  # errors
  expect_error(expand_hesim_data(hesim.dat, by = c("strategies", "patients", 
                                                  "states", "transitions")))
  expect_error(expand_hesim_data(hesim.dat, by = c("strategies", "patients", 
                                                  "states", "wrong_table")))
  hesim.dat2 <- hesim.dat[c("strategies", "patients")]
  expect_error(expand_hesim_data(hesim.dat2, by = c("strategies", "patients", 
                                                  "states")))
})

# input_data class -------------------------------------------------------------
# By treatment strategy and patient
dat <- expand_hesim_data(hesim.dat)$data
input.dat <- input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)))

## Number of rows in X is inconsistent with strategy_id 
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id[-1],
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id))))

## Size of patient_id is incorrect
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = sort(dat$patient_id),
                       n_patients = length(unique(dat$strategy_id))))

## n_patients is incorrect v1
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$strategy_id,
                       n_patients = length(unique(dat$strategy_id))))

## n_patients is incorrect v2
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = 1))

## patient_id is not sorted correctly
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = sort(dat$patient_id),
                       n_patients = length(unique(dat$patient_id))))


# By treatment strategy, line, and patient
dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", "lines"))$data
n.lines <- hesim.dat$lines[, .N, by = "strategy_id"]
input.dat <- input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)),
                       line = dat$line,
                       n_lines = n.lines)

## n_lines is incorrect v1
n.lines[, N := N + 1]
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)),
                       line = dat$line,
                       n_lines = n.lines))

## n_lines is incorrect v2
n.lines[, N := N - 1]
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)),
                       line = dat$line,
                       n_lines = n.lines + 1))

## line is not sorted correctly
expect_error(input_data(X = model.matrix(~ age, dat),
                       strategy_id = dat$strategy_id,
                       n_strategies = length(unique(dat$strategy_id)),
                       patient_id = dat$patient_id,
                       n_patients = length(unique(dat$patient_id)),
                       line = sort(dat$line),
                       n_lines = n.lines))

# form_input_data with formula objects -----------------------------------------
test_that("form_input_data.formula", {
  f <- formula(~ age)
  dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", "lines", "states"))
  input.dat <- form_input_data(f, dat)
  
  expect_equal(input.dat$state_id, dat$data$state_id)
  expect_equal(input.dat$strategy_id, dat$data$strategy_id)
  expect_equal(input.dat$line, dat$data$line)
  expect_equal(input.dat$patient_id, dat$data$patient_id)
  expect_equal(ncol(input.dat$X), 2)
  expect_equal(as.numeric(input.dat$X[, "age"]), dat$data$age)
  
  class(dat) <- "data.table"
  expect_error(form_input_data(f, dat, 
                              id_vars = c("strategy_id", "patient_id")))
})

test_that("form_input_data.formula_list", {
  dat <- expand_hesim_data(hesim.dat)
  f.list <- formula_list(list(f1 = formula(~ age), f2 = formula(~ 1)))
  expect_equal(class(f.list), "formula_list")
  input.dat <- form_input_data(f.list, dat)
  
  expect_equal(length(input.dat$X), length(f.list))
  expect_equal(names(input.dat$X), names(f.list))
  expect_equal(as.numeric(input.dat$X$f1[, "age"]), dat$data$age)
  expect_equal(ncol(input.dat$X$f1), 2)
  expect_equal(ncol(input.dat$X$f2), 1)
})

# form_input_data with lm objects ----------------------------------------------
dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", "states"))
fit1 <- stats::lm(costs ~ female + state_name, data = part_surv4_simdata$costs$medical)

test_that("form_input_data.lm", {
  input.dat <- form_input_data(fit1, dat)
  
  expect_equal(ncol(input.dat$X), 4)
  expect_equal(as.numeric(input.dat$X[, "female"]), dat$data$female)
})

test_that("form_input_data.lm_list", {
  fit2 <- stats::lm(costs ~ 1, data = part_surv4_simdata$costs$medical)
  fit.list <- hesim:::lm_list(fit1 = fit1, fit2 = fit2)
  input.dat <- form_input_data(fit.list, dat)
  
  expect_equal(ncol(input.dat$X$fit1), 4)
  expect_equal(ncol(input.dat$X$fit2), 1)
  expect_equal(as.numeric(input.dat$X$fit1[, "female"]), dat$data$female)
})

# form_input_data with flexsurvreg objects -------------------------------------
test_that("form_input_data.flexsurv", {
  dat <- expand_hesim_data(hesim.dat)
  fit <- flexsurv::flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                     anc = list(sigma = ~ group), dist = "gengamma") 
  input.dat <- form_input_data(fit, dat)
  
  expect_equal(input.dat$strategy_id, dat$data$strategy_id)
  expect_equal(input.dat$state_id, dat$data$state_id)
  expect_equal(input.dat$patient_id, dat$data$patient_id)
  expect_equal(class(input.dat$X), "list")
  expect_equal(class(input.dat$X[[1]]), "matrix")
  expect_equal(length(input.dat$X), 3)
  expect_equal(ncol(input.dat$X$mu), 3)
  expect_equal(ncol(input.dat$X$sigma), 3)
  expect_equal(ncol(input.dat$X$Q), 1)
})

fit1.wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                              data = ovarian, dist = "weibull")
fit1.exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                              data = ovarian, dist = "exp")
flexsurvreg.list1 <- flexsurvreg_list(wei = fit1.wei, exp = fit1.exp)
dat <- expand_hesim_data(hesim.dat)

test_that("form_input_data.flexsurv_list", {
  input.dat <- form_input_data(flexsurvreg.list1, dat)  
  
  expect_equal(class(input.dat$X$wei$shape), "matrix")
})

fit2.wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1 + age, 
                              data = ovarian, dist = "weibull")
fit2.exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1 + age, 
                              data = ovarian, dist = "exp")
flexsurvreg.list2 <- flexsurvreg_list(wei = fit2.wei, exp = fit2.exp)
joined.flexsurvreg.list <- joined_flexsurvreg_list(mod1 = flexsurvreg.list1,
                                                   mod2 = flexsurvreg.list2,
                                                   times = list(2, 5))

test_that("form_input_data.joined_flexsurv_list", {
  input.dat <- form_input_data(joined.flexsurvreg.list, dat)  
  
  expect_equal(input.dat$state_id, dat$state_id)
  expect_equal(class(input.dat$X[[1]]$wei$shape), "matrix")
})
