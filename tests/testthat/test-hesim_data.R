context("input_data unit tests")
library("flexsurv")
library("data.table")
rm(list = ls())

strategies <- data.table(
  strategy_id = c(1, 2),
  strategy_name = c("Strategy 1", "Strategy 2")
)
patients <- data.table(
  patient_id = 1:3, 
  age = c(45, 47, 60),
  female = c(1, 0, 0),
  grp_id = 1:3,
  group = factor(c("Good", "Medium", "Poor"))
) 
states <- data.frame(
  state_id =  seq(1, 3),
  state_name = factor(paste0("state", seq(1, 3)))
)
trans <- data.frame(
  transition_id = seq(1, 4),
  from = c(1, 1, 2, 2),
  to = c(2, 3, 1, 3),
  transition_name = c("1->2", "1->3", "2->1", "2->3")
)
hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients,
  states = states,
  transitions = trans
)

# create_lines_dt --------------------------------------------------------------
test_that("create_lines_dt() works as expected", {
  lines_dt <- create_lines_dt(list(c(1, 2, 5), c(1, 2)))
  
  expect_true(inherits(lines_dt, "data.table"))
  expect_equal(lines_dt$treatment_id[3], 5)
  expect_equal(lines_dt$line, 
               c(seq(1, 3), seq(1, 2)))
  
  # explicit strategy ids
  lines_dt <- create_lines_dt(list(c(1, 2, 5), c(1, 2)),
                              strategy_ids = c(3, 5))
  expect_equal(lines_dt$strategy_id, c(3, 3, 3, 5, 5))
})

test_that("create_lines_dt() throws error if the 'strategy_list' does not contain integers ", {
  expect_error(
    create_lines_dt(list(c("tx1", "tx2"), c("tx1"))),
    "Elements in 'strategy_list' should be integers."
  )
})


# create_trans_dt() ------------------------------------------------------------
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))

test_that("create_trans_dt() returns correct columns and values", {
  trans_dt <- create_trans_dt(tmat)
  
  expect_true(inherits(trans_dt, "data.table"))
  expect_equal(trans_dt$transition_id, 
               c(1, 2, 3))
  expect_equal(trans_dt$from, 
               c(1, 1, 2))
  expect_equal(trans_dt$to, 
               c(2, 3, 3))
})
  
test_that(paste0("create_trans_dt() does not create '_name' columns without both ",
                 "rownames and column names for 'trans_mat'"), {
  rownames(tmat) <- c("No BOS", "BOS", "Death")
  trans_dt <- create_trans_dt(tmat)
  expect_equal(trans_dt$from_name, NULL)
})

test_that("create_trans_dt() automatically creates '_name' columns", {
  colnames(tmat) <- rownames(tmat) <- c("No BOS", "BOS", "Death")
  trans_dt <- create_trans_dt(tmat)
  expect_equal(trans_dt$from_name, rownames(tmat)[c(1, 1, 2)])
  expect_equal(trans_dt$to_name, colnames(tmat)[c(2, 3, 3)])
})

# hesim data() -----------------------------------------------------------------
test_that("hesim_data() works with strategies and patients", {
  h <- hesim_data(strategies = strategies,
                  patients = patients)
  
  expect_true(inherits(h, "hesim_data"))
  expect_equal(h$state, NULL)
  expect_equal(h$patients, patients)
})

test_that("hesim_data() works with strategies, patients, and states", {
  h <- hesim_data(strategies = strategies,
                  patients = patients, 
                  states = states)
  expect_equal(h$states, states)
  
})

test_that("expand.hesim_data() works with strategies", {
  e <- expand(hesim_dat, by = c("strategies"))
  expect_equal(e, data.table(strategies), check.attributes = FALSE)
  expect_equal(attributes(e)$id_vars, "strategy_id")
  
})

test_that("expand.hesim_data() works with strategies and patients", {
  e <- expand(hesim_dat, by = c("strategies", "patients"))
  expect_equal(nrow(e), 
               nrow(strategies) * nrow(patients))
  expect_equal(attributes(e)$id_vars, c("strategy_id", "patient_id"))
})

test_that("Order of 'by' in expand.hesim_data() does not matter", {
  e1 <- expand(hesim_dat, by = c("strategies", "patients"))
  e2 <- expand(hesim_dat, by = c("strategies", "patients"))
  expect_equal(e1, e2)
})

test_that("expand.hesim_data() works with strategies, patients, and time intervals", {
  e <- expand(hesim_dat, by = c("strategies", "patients"),
              times = c(0, 2, 4))
  expect_equal(nrow(e), 
               nrow(strategies) * nrow(patients) * 3)
})

test_that("expand.hesim_data() throws error if both 'states' and 'transitions; are in 'by'", {
  expect_error(
    expand(hesim_dat, by = c("strategies", "patients", "states", "transitions")),
    "Cannot expand by both 'transitions' and 'states'."
  )
})

test_that("expand.hesim_data() throws error if incorrect table is in 'by'", {
  expect_error(
    expand(hesim_dat, by = c("strategies", "patients", "states", "wrong_table")),
    "One of the elements specified in 'by' is not a table in 'hesim_data'."
  )
})

test_that("expand.hesim_data() throws error if element table is in 'by' is not in hesim data", {
  h <- hesim_dat[c("strategies", "patients")]
  class(h) <-"hesim_data"
  expect_error(
    expand(h, by = c("strategies", "patients", "states")),
    "Cannot merge a NULL data table."
  )
})

test_that("expand.hesim_data() preserves attributes when subsetting", {
  
  # with data table
  dat <- expand(hesim_dat)
  expect_equal(attributes(dat[1])$id_vars, c("strategy_id", "patient_id"))
  expect_equal(dat[1:2, age], hesim_dat$patients$age[1:2], check.attributes = FALSE)
  tmp <- dat[1:2, .(age, female)]
  expect_equal(nrow(tmp), 2)
  expect_equal(colnames(tmp), c("age", "female"))
  expect_equal(attributes(tmp)$id_vars, c("strategy_id", "patient_id"))
  
  # with data frame
  setattr(dat, "class", c("expanded_hesim_data", "data.frame"))
  expect_equal(attributes(dat[1, ])$id_vars, c("strategy_id", "patient_id"))
  tmp <- dat[, c("age", "female")]
  expect_equal(nrow(tmp), nrow(dat))
  expect_equal(colnames(tmp), c("age", "female"))
  expect_equal(attributes(tmp)$id_vars, c("strategy_id", "patient_id"))
})

# ID attributes ----------------------------------------------------------------
d1 <- expand(hesim_dat)
d2 <- expand(hesim_dat, by = c("strategies", "patients", "states"))

test_that("id_attributes() works", {
  # Treatment strategies and patients
  id <- id_attributes(strategy_id = d1$strategy_id,
                      n_strategies = length(unique(d1$strategy_id)),
                      patient_id = d1$patient_id,
                      n_patients = length(unique(d1$patient_id)))
  expect_true(inherits(id, "id_attributes"))
  
  # Treatment strategies, patients, and health state
  id <- id_attributes(strategy_id = d2$strategy_id,
                      n_strategies = length(unique(d2$strategy_id)),
                      patient_id = d2$patient_id,
                      n_patients = length(unique(d2$patient_id)),
                      state_id = d2$state_id,
                      n_states = length(unique(d2$state_id)))
  expect_true(inherits(id, "id_attributes"))
})

test_that("id_attributes() throws error if the size of an attribute is wrong ", {
  expect_error(
    id_attributes(strategy_id = d1$strategy_id,
                  n_strategies = length(unique(d1$strategy_id)),
                  patient_id = d1$patient_id,
                  n_patients = length(unique(d1$patient_id)) + 2),
    "The number of unique observations in 'patient_id' does not equal 'n_patients'."
  )
})

test_that("id_attributes() throws error if length of ID vectors is wrong ", {
  expect_error(
    id_attributes(strategy_id = d2$strategy_id,
                  n_strategies = length(unique(d2$strategy_id)),
                  patient_id = d2$patient_id,
                  n_patients = length(unique(d2$patient_id))),
    paste0("The length of the ID variables is not consistent with the number of ",
           "unique values of each ID variable.")
  )
})

test_that("id_attributes() throws error if the sorting order is wrong ", {
  # Treatment strategies and patients
  d <- copy(d1)
  setorderv(d, c("patient_id", "strategy_id"))
  expect_error(
    id_attributes(strategy_id = d$strategy_id,
                  n_strategies = length(unique(d$strategy_id)),
                  patient_id = d$patient_id,
                  n_patients = length(unique(d$patient_id))),
    paste0("The ID variables are not sorted correctly. The sort priority of the ", 
            "ID variables must be as follows: strategy_id and patient_id.")
  )
  
  # Treatment strategies, patients, and health state
  d <- copy(d2)
  setorderv(d, c("strategy_id", "state_id", "patient_id"))
  expect_error(
    id_attributes(strategy_id = d$strategy_id,
                  n_strategies = length(unique(d$strategy_id)),
                  patient_id = d$patient_id,
                  n_patients = length(unique(d$patient_id)),
                  state_id = d$state_id,
                  n_states = length(unique(d$state_id))),
    paste0("The ID variables are not sorted correctly. The sort priority of the ", 
           "ID variables must be as follows: strategy_id, patient_id, and state_id.")
  )
})

test_that("id_attributes() throws error if there are not the right sizes within groups ", {
  # Treatment strategies and patients
  d <- copy(d1)
  d[, patient_id := ifelse(strategy_id == 1 & patient_id == 2, 
                         1, patient_id)]
  expect_error(
    id_attributes(strategy_id = d$strategy_id,
                  n_strategies = length(unique(d$strategy_id)),
                  patient_id = d$patient_id,
                  n_patients = length(unique(d$patient_id))),
    paste0("The number of unique patient_id observations within each strategy_id ",
           "group must equal n_patients.")
  )
  
  # Treatment strategies, patients, and health state
  d <- copy(d2)
  d[, state_id := ifelse(strategy_id == 1 & patient_id == 1 & state_id == 2, 
                         1, state_id)]
  expect_error(
    id_attributes(strategy_id = d$strategy_id,
                  n_strategies = length(unique(d$strategy_id)),
                  patient_id = d$patient_id,
                  n_patients = length(unique(d$patient_id)),
                  state_id = d$state_id,
                  n_states = length(unique(d$state_id))),
    paste0("The number of unique state_id observations within each strategy_id ",
           "and patient_id group must equal n_states.")
  )
})

# get_labels() -----------------------------------------------------------------
test_that("get_labels() works with 4 tables in hesim_data", {
  x <- get_labels(hesim_dat, grp = "group", death_label = NULL)
  expect_true(is.list(x))
  expect_equal(length(x), length(hesim_dat))
  expect_equivalent(sapply(x, length),
                    sapply(hesim_dat, nrow))
})

test_that("get_labels() adds a death state if death_label is not NULL", {
  x <- get_labels(hesim_dat, grp = "group")
  expect_true("Death" %in% names(x$state_id))
})

test_that("get_labels() works with 2 tables in hesim_data", {
  h <- hesim_data(strategies = strategies, patient = patients)
  x <- get_labels(h, grp = "group")
  
  expect_equivalent(x[[1]], h$strategies$strategy_id)
  expect_equivalent(names(x[[1]]), h$strategies$strategy_name)
  
  expect_equivalent(x[[2]], h$patients$grp_id)
  expect_equivalent(names(x[[2]]), as.character(h$patients$group))
})

test_that("get_labels() works with only 1 label", {
  x <- get_labels(hesim_dat,  grp = NULL, state = NULL, transition = NULL)
  expect_equal(length(x), 1)
})

test_that("get_labels() works with more patients than subgroups", {
  pt <- data.table(patient_id = 1:3, grp_id = c(1, 2, 2), grp_name = c("g1", "g2", "g2"))
  h <- hesim_data(strategies = strategies, patient = pt)
  x <- get_labels(h)
  expect_equivalent(x$grp_id, unique(pt$grp_id))
  expect_equivalent(names(x$grp_id), as.character(unique(pt$grp_name)))
})

test_that("get_labels() removes label if variable does not exist", {
  x <- get_labels(hesim_dat)
  expect_equal(names(x), c("strategy_id", "state_id", "transition_id"))
})


test_that("get_labels() throws an error if there is not exactly one label for each ID", {
  pt <- data.table(patient_id = 1:3, grp_id = c(1, 2, 2), grp_name = c("g1", "g2", "g3"))
  h <- hesim_data(strategies = strategies, patient = pt)
  expect_error(
    get_labels(h), 
    "There should be exactly one label for each ID value."
  )
})

test_that("get_labels() throws error with all labels NULL", {
  expect_error(
    get_labels(hesim_dat, strategy = NULL, grp = NULL, 
              state = NULL, transition = NULL),
    "There are no labels to get."
  )
})

test_that("get_labels() throws error with no valid labels", {
  expect_error(
    get_labels(hesim_dat, strategy = "s", grp = "g",
               state = "s2", transition ="t"),
    "The selected labels are not contained in the tables of 'object'."
  )
})

# set_labels() -----------------------------------------------------------------
labs <- get_labels(hesim_dat, grp = "group")
d <- data.table(strategy_id = rep(1:2, each = 3) , grp_id = rep(1:3, 2))

test_that("set_labels() modifies existing variables if 'new_names' = NULL", {
  d2 <- copy(d)
  set_labels(d2, labels = labs, as_factor = FALSE)
  expect_equal(unique(d2$strategy_id), names(labs$strategy_id))
  expect_equal(unique(d2$grp_id), names(labs$grp_id))
})

test_that("set_labels() creates new variables if 'new_names' = NULL", {
  d2 <- copy(d)
  set_labels(d2, labels = labs, new_names = c("s", "g"), as_factor = FALSE)
  expect_equal(unique(d2$s), names(labs$strategy_id))
  expect_equal(unique(d2$g), names(labs$grp_id))
})

test_that("set_labels() creates factors if 'as_factor' = TRUE", {
  d2 <- copy(d)
  set_labels(d2, labels = labs, as_factor = TRUE)
  expect_equal(levels(d2$strategy_id), names(labs$strategy_id))
  expect_equal(levels(d2$grp_id), names(labs$grp_id))
})

test_that("set_labels() does nothing if there are no labels", {
  d2 <- copy(d)
  set_labels(d2, labels = NULL)
  expect_equal(d, d2)
})