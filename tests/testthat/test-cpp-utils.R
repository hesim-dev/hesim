context("util.cpp unit tests")

# Test c++ function add_constant -----------------------------------------------
test_that("add_constant", {
  v1 <- c(1, 2, 3)
  v2 <- c(1.2, 4, 5.6)
  expect_equal(hesim:::C_test_add_constant_int(v1, 5.5),
               floor(v1 + 5.5))
  expect_equal(hesim:::C_test_add_constant_double(v2, 5.5),
               v2 + 5.5)
})

# Test c++ function pv ---------------------------------------------------------
test_that("pv", {
  expect_equal(hesim:::C_test_pv(1, 0, 0, 4),
               4)
  expect_equal(hesim:::C_test_pv(2, .03, 1, 4),
               2 * ((exp(-.03 * 1) - exp(-.03 * 4))/.03))
})

# Test c++ function seq --------------------------------------------------------
test_that("seq", {
  expect_equal(hesim:::C_test_seq(0, 3.01, .1),
               seq(0, 3, .1))
  expect_equal(hesim:::C_test_seq(0, 2.9, .2),
               seq(0, 2.9, .2))
  expect_equal(hesim:::C_test_seq(0, -3.01, -.1),
               seq(0, -3, -.1))
  expect_equal(hesim:::C_test_seq(2, -2.45, -.13),
               seq(2, -2.45, -.13))
  expect_error(hesim:::C_test_seq(0, -3, .1))
  expect_error(hesim:::C_test_seq(0, 3, -.1))
})

# Test c++ function max_lt -----------------------------------------------------
test_that("max_lt", {
  expect_equal(hesim:::C_test_max_lt(seq(0, 10), 5), 5)
  expect_equal(hesim:::C_test_max_lt(seq(0, 10), 5.1), 5)
  expect_equal(hesim:::C_test_max_lt(seq(0, 10), Inf), 10)
})

# Test c++ function hesim_bound -----------------------------------------------------
test_that("find_interval", {
    vec = c(0,3,5)
    expect_equal(hesim:::C_test_find_interval(-1,vec),0)
    expect_equal(hesim:::C_test_find_interval(0,vec),0)
    expect_equal(hesim:::C_test_find_interval(1,vec),0)
    expect_equal(hesim:::C_test_find_interval(3,vec),1)
    expect_equal(hesim:::C_test_find_interval(3.1,vec),1)
    expect_equal(hesim:::C_test_find_interval(5,vec),2)
    expect_equal(hesim:::C_test_find_interval(5.1,vec),2)
})

# Test c++ function member_of -----------------------------------------------------
test_that("member_of", {
    lookup = c(0,3,5)
    expect_equal(hesim:::C_test_is_member_of(-1,lookup),FALSE)
    expect_equal(hesim:::C_test_is_member_of(0,lookup),TRUE)
    expect_equal(hesim:::C_test_is_member_of(0.1,lookup),FALSE)
    expect_equal(hesim:::C_test_is_member_of(3,lookup),TRUE)
    expect_equal(hesim:::C_test_is_member_of(5,lookup),TRUE)
    expect_equal(hesim:::C_test_is_member_of(6,lookup),FALSE)
})
