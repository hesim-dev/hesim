context("utils.R unit tests")

# list depth -------------------------------------------------------------------
list1 <- list(1)
list2 <- list(1, list1)
list3 <- list(1, list1, list2)

expect_equal(list_depth(list1), 1)
expect_equal(list_depth(list2), 2)
expect_equal(list_depth(list3), 3)
