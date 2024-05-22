context("inline unit tests")
library("Rcpp")

# Partitioned survival fits  ---------------------------------------------------
test_that("Show that inline works for a simple example ", {
    sourceCpp(code="
// [[Rcpp::depends(hesim)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <hesim.h>
// [[Rcpp::export]]
double test_inline_gengamma(double mu, double sigma, double Q) {
  hesim::stats::gengamma gg(mu, sigma, Q);
  return gg.random();
}")
    set.seed(12345)
    expect_true(abs(test_inline_gengamma(1.0, 1.0, 1.0) - 2.717582) < 1e-5)
})
