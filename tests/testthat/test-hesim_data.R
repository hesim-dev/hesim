context("input_data unit tests")
library("flexsurv")
library("data.table")
rm(list = ls())

strategies_dt <- data.table(
  strategy_id = c(1, 2),
  strategy_name = c("Strategy 1", "Strategy 2")
)
patients_dt <- data.table(
  patient_id = 1:3, 
  age = c(45, 47, 60),
  female = c(1, 0, 0),
  grp_id = 1:3,
  group = factor(c("Good", "Medium", "Poor"))
) 
states_dt <- data.frame(
  state_id =  seq(1, 3),
  state_name = factor(paste0("state", seq(1, 3)))
)
trans_dt <- data.frame(
  transition_id = seq(1, 4),
  from = c(1, 1, 2, 2),
  to = c(2, 3, 1, 3),
  transition_name = c("1->2", "1->3", "2->1", "2->3")
)
hesim_dat <- hesim_data(
  strategies = strategies_dt,
  patients = patients_dt,
  states = states_dt,
  transitions = trans_dt
)

# create_lines_dt --------------------------------------------------------------
test_that("create_lines_dt", {
  lines_dt <- create_lines_dt(list(c(1, 2, 5), c(1, 2)))
  
  expect_true(inherits(lines_dt, "data.table"))
  expect_equal(lines_dt$treatment_id[3], 5)
  expect_equal(lines_dt$line, 
               c(seq(1, 3), seq(1, 2)))
  
  # explicit strategy ids
  lines_dt <- create_lines_dt(list(c(1, 2, 5), c(1, 2)),
                              strategy_ids = c(3, 5))
  expect_equal(lines_dt$strategy_id, c(3, 3, 3, 5, 5))
  
  # errors
  expect_error(create_lines_dt(list(c("tx1", "tx2"),
                                  c("tx1"))))
})

# create_trans_dt() ------------------------------------------------------------
test_that("create_trans_dt", {
  tmat <- rbind(c(NA, 1, 2),
                c(NA, NA, 3),
                c(NA, NA, NA))
  trans_dt <- create_trans_dt(tmat)
  
  expect_true(inherits(trans_dt, "data.table"))
  expect_equal(trans_dt$transition_id, 
               c(1, 2, 3))
  expect_equal(trans_dt$from, 
               c(1, 1, 2))
  expect_equal(trans_dt$to, 
               c(2, 3, 3))
  
  # Row and column names
  rownames(tmat) <- c("No BOS", "BOS", "Death")
  trans_dt <- create_trans_dt(tmat)
  expect_equal(trans_dt$from_name, NULL)
  
  colnames(tmat) <- rownames(tmat)
  trans_dt <- create_trans_dt(tmat)
  expect_equal(trans_dt$from_name, rownames(tmat)[c(1, 1, 2)])
  expect_equal(trans_dt$to_name, colnames(tmat)[c(2, 3, 3)])
})

# hesim data() -----------------------------------------------------------------
test_that("hesim_data", {

  # strategy by patient
  hesim_dat2 <- hesim_data(strategies = strategies_dt,
                          patients = patients_dt)
  
  expect_true(inherits(hesim_dat2, "hesim_data"))
  expect_equal(hesim_dat2$state, NULL)
  expect_equal(hesim_dat2$patients, patients_dt)
  
  # strategy by patient by state
  hesim_dat2 <- hesim_data(strategies = strategies_dt,
                            patients = patients_dt, 
                            states = states_dt)
  expect_equal(hesim_dat2$states, states_dt)
  
  # Expand
  expanded_dt <- expand(hesim_dat, by = c("strategies"))
  expect_equal(expanded_dt, data.table(strategies_dt), check.attributes = FALSE)
  expect_equal(attributes(expanded_dt)$id_vars, "strategy_id")
  
  expanded_dt <- expand(hesim_dat, by = c("strategies", "patients"))
  expanded_dt2 <- expand(hesim_dat, by = c("patients", "strategies"))
  expect_equal(nrow(expanded_dt), 
               nrow(strategies_dt) * nrow(patients_dt))
  expect_equal(expanded_dt, expanded_dt2)
  expect_equal(attributes(expanded_dt)$id_vars, attributes(expanded_dt2)$id_vars)
  expect_equal(attributes(expanded_dt)$id_vars, c("strategy_id", "patient_id"))
  
  expanded_dt <- expand(hesim_dat, by = c("strategies", "patients"),
                        times = c(0, 2, 4))
  expect_equal(nrow(expanded_dt), 
               nrow(strategies_dt) * nrow(patients_dt) * 3)
  
  
  # errors
  expect_error(expand(hesim_dat, by = c("strategies", "patients", 
                                                  "states", "transitions")))
  expect_error(expand(hesim_dat, by = c("strategies", "patients", 
                                                  "states", "wrong_table")))
  hesim_dat2 <- hesim_dat[c("strategies", "patients")]
  class(hesim_dat2) <-"hesim_data"
  expect_error(expand(hesim_dat2, by = c("strategies", "patients", 
                                                  "states")))
  
  # Attributes are preserved with subsetting
  ## with data table
  dat <- expand(hesim_dat)
  expect_equal(attributes(dat[1])$id_vars, c("strategy_id", "patient_id"))
  expect_equal(dat[1:2, age], hesim_dat$patients$age[1:2], check.attributes = FALSE)
  tmp <- dat[1:2, .(age, female)]
  expect_equal(nrow(tmp), 2)
  expect_equal(colnames(tmp), c("age", "female"))
  expect_equal(attributes(tmp)$id_vars, c("strategy_id", "patient_id"))
  
  ## with data frame
  setattr(dat, "class", c("expanded_hesim_data", "data.frame"))
  expect_equal(attributes(dat[1, ])$id_vars, c("strategy_id", "patient_id"))
  tmp <- dat[, c("age", "female")]
  expect_equal(nrow(tmp), nrow(dat))
  expect_equal(colnames(tmp), c("age", "female"))
  expect_equal(attributes(tmp)$id_vars, c("strategy_id", "patient_id"))
})

# ID attributes ----------------------------------------------------------------
# By treatment strategy and patient
dat <- expand(hesim_dat)
id <- id_attributes(strategy_id = dat$strategy_id,
                    n_strategies = length(unique(dat$strategy_id)),
                    patient_id = dat$patient_id,
                    n_patients = length(unique(dat$patient_id)))
expect_true(class(id) == "id_attributes")

## Size of patient_id is incorrect
expect_error(id_attributes(strategy_id = dat$strategy_id,
                           n_strategies = length(unique(dat$strategy_id)),
                           patient_id = dat$patient_id,
                           n_patients = length(unique(dat$strategy_id))))

## n_patients is incorrect 
expect_error(id_attributes(strategy_id = dat$strategy_id,
                           n_strategies = length(unique(dat$strategy_id)),
                           patient_id = dat$patient_id,
                           n_patients = 1))

## patient_id is not sorted correctly
expect_error(id_attributes(strategy_id = dat$strategy_id,
                           n_strategies = length(unique(dat$strategy_id)),
                           patient_id = sort(dat$patient_id),
                           n_patients = length(unique(dat$patient_id))))

# get_labels() -----------------------------------------------------------------
test_that("get_labels() works with 4 tables in hesim_data", {
  x <- get_labels(hesim_dat, grp_label = "group")
  expect_true(is.list(x))
  expect_equal(length(x), length(hesim_dat))
  expect_equivalent(sapply(x, length),
                    sapply(hesim_dat, nrow))
})

test_that("get_labels() works with 2 tables in hesim_data", {
  h <- hesim_data(strategies = strategies_dt, patient = patients_dt)
  x <- get_labels(h, grp_label = "group")
  
  expect_equivalent(x[[1]], h$strategies$strategy_id)
  expect_equivalent(names(x[[1]]), h$strategies$strategy_name)
  
  expect_equivalent(x[[2]], h$patients$grp_id)
  expect_equivalent(names(x[[2]]), as.character(h$patients$group))
})

test_that("get_labels() works with only 1 label", {
  x <- get_labels(hesim_dat,  grp_label = NULL, state_label = NULL, transition_label = NULL)
  expect_equal(length(x), 1)
})

test_that("get_labels() works with more patients than subgroups", {
  pt <- data.table(patient_id = 1:3, grp_id = c(1, 2, 2), grp_name = c("g1", "g2", "g2"))
  h <- hesim_data(strategies = strategies_dt, patient = pt)
  x <- get_labels(h)
  expect_equivalent(x$grp_id, unique(pt$grp_id))
  expect_equivalent(names(x$grp_id), as.character(unique(pt$grp_name)))
})

test_that("get_labels() throws an error if there is not exactly one label for each ID", {
  pt <- data.table(patient_id = 1:3, grp_id = c(1, 2, 2), grp_name = c("g1", "g2", "g3"))
  h <- hesim_data(strategies = strategies_dt, patient = pt)
  expect_error(
    get_labels(h), 
    "There should be exactly one label for each ID value."
  )
})

test_that("get_labels() throws error with no labels", {
  expect_error(
    get_labels(hesim_dat, strategy_label = NULL,  grp_label = NULL, 
              state_label = NULL, transition_label = NULL),
    "There are no labels to get."
  )
})

test_that("get_labels() throws error if label variable does not exist", {
  expect_error(get_labels(hesim_dat),
               "'grp_name' is not contained in the 'patients' table")
})

