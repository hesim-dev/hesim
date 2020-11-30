context("data.R unit tests")
library("data.table")

# Convert multi-state dataset to PFS/OS dataset --------------------------------
msdata <- data.table(
  patient_id = c(rep(1, 2),
                 rep(2, 3)),
  transition_id = c(1, 2, 
                    1, 2, 3),
  status = c(0, 1,
             1, 0, 0),
  tstop = c(10, 10,
            5, 5, 15)
)
msdata_pfs_os <- as_pfs_os(msdata, patient_vars = "patient_id", time_stop = "tstop")

test_that("as_pfs_os() returns correct output given stable -> death trajectory" , {
  y <- msdata_pfs_os[patient_id == 1]
  expect_equal(y$pfs_time, 10)
  expect_equal(y$pfs_status, 1)
  expect_equal(y$os_time, 10)
  expect_equal(y$os_status, 1)
})

test_that(paste0("as_pfs_os() returns correct output given stable -> progression -> death trajectory",
                 " and right censoring for progression -> death transition"), {
  y <- msdata_pfs_os[patient_id == 2]
  expect_equal(y$pfs_time, 5)
  expect_equal(y$pfs_status, 1)
  expect_equal(y$os_time, 15)
  expect_equal(y$os_status, 0)
})

test_that(paste0("as_pfs_os() returns correct output given stable -> progression -> death trajectory",
                 " and observed progression -> death transition"), {
  x <- data.frame(
    ptid = rep(1, 3),
    trans = c(1, 2, 3),
    status = c(1, 0, 1),
    time_stop = c(7, 7, 11)
   )
   y <- as_pfs_os(x, patient_vars = "ptid", transition = "trans")
   expect_equal(y$pfs_time, 7)
   expect_equal(y$pfs_status, 1)
   expect_equal(y$os_time, 11)
   expect_equal(y$os_status, 1)
})

test_that("as_pfs_os() returns correct output given stable -> death is right censored" , {
  x <- data.frame(
    ptid = rep(1, 2),
    transition_id = c(1, 2),
    died = c(0, 0),
    time_stop = c(5, 5)
  )
  y <- as_pfs_os(x, patient_vars = "ptid", status = "died")
  expect_equal(y$pfs_time, 5)
  expect_equal(y$os_time, 5)
  expect_equal(y$pfs_status, 0)
  expect_equal(y$os_status, 0)
})

test_that("as_pfs_os() returns correct number of rows" , {
  expect_equal(
    nrow(msdata_pfs_os), 
    length(unique(msdata$patient_id))
  )
})

test_that("as_pfs_os() returns error with wrong transition variable" , {
  expect_error(
    as_pfs_os(onc3, patient_vars = "age", transition = "patient_id"),
    "'patient_id' should be a vector with unique values c(1, 2, 3).",
    fixed = TRUE
    )
})

test_that("as_pfs_os() returns error with transition numbers not in c(1, 2, 3)" , {
  y <- copy(msdata)[, transition_id := ifelse(transition_id == 3, 4, transition_id)]
  expect_error(
    as_pfs_os(y, patient_vars = "patient_id", time_stop = "tstop"),
    "'transition_id' should be a vector with unique values c(1, 2, 3).",
    fixed = TRUE
  )
})

test_that("as_pfs_os() returns error with only one transition", {
  expect_error(
    as_pfs_os(msdata[transition_id == 1], patient_vars = "patient_id", 
              time_stop = "tstop"),
    "'transition_id' should contain at a minimum values 1 and 2.",
  )
})

test_that("as_pfs_os() returns error with only two transition that aren't 1 and 2", {
  expect_error(
    as_pfs_os(msdata[transition_id %in% c(1, 3)], patient_vars = "patient_id", 
              time_stop = "tstop"),
    "If 'transition_id' contains 2 values, they should be 1 and 2.",
  )
})