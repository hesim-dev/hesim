context("quad.h unit tests")

expect_equal(hesim:::test_quad_functor(), 2 * 2)

expect_equal(hesim:::test_quad_lambda(), 1/3 *  2^3)

