test_that("`.var()` returns expected results", {
  
  expect_equal(.var(1.0), 0.0)
  expect_equal(.var(rep(1.0, 10L)), 0.0)
  
  expected <- 2.0 * (4.5^2 + 3.5^2 + 2.5^2 + 1.5^2 + 0.5^2) / 10
  expect_equal(.var(1:10), expected)
  
  expect_equal(.var(1:10), var(1:10) * 9.0 / 10.0)
  
  expected <- 0.0
  for (i in 1:10) expected <- expected + i^2
  expect_equal(.var(1:10, 5L), expected / 5)

})

test_that("`.calEq()` returns expected results, intercept only", {
  
  X <- matrix(0.0, nrow = 100L, ncol = 0L)
  par <- 0.25
  
  expect_equal(.calEqs(par, X, X), exp(0.25)*100.0 - 100)
  
})

test_that("`.calEq()` returns expected results, 1 covariate", {

  X <- withr::with_seed(1234L, matrix(stats::rnorm(100), ncol = 1L))
  par <- c(0.25, 0.5)
  wgt <- drop(exp(0.25 + X * 0.5))
  
  expect_equal(.calEqs(par, X, X), 
               c(sum(wgt), sum(wgt * X)) - c(100, sum(X)))
  
})

test_that("`.calEq()` returns expected results, >1 covariate", {
  
  X <- withr::with_seed(1234L, matrix(stats::rnorm(100), ncol = 2L))
  par <- c(0.25, 0.5, 0.75)
  wgt <- drop(exp(0.25 + 0.5 * X[, 1L] + 0.75 * X[, 2L]))
  
  expect_equal(.calEqs(par, X, X), 
               c(sum(wgt), sum(wgt * X[, 1L]), sum(wgt * X[, 2L])) - 
                 c(50, sum(X[, 1L]), sum(X[, 2L])))
  
})

test_that("`.calEq()` returns expected results, >1 covariate", {
  
  X <- withr::with_seed(1234L, matrix(stats::rnorm(100), ncol = 2L))
  par <- c(0.25, 0.5, 0.75)
  wgt <- drop(exp(0.25 + 0.5 * X[, 1L] + 0.75 * X[, 2L]))
  
  expect_equal(.calEqs(par, 2.0*X, X), 
               c(sum(wgt), sum(wgt * X[, 1L]), sum(wgt * X[, 2L])) - 
                 c(50, 2.0*sum(X[, 1L]), 2.0*sum(X[, 2L])))
  
})

test_that("`.calEqGMM()` returns expected results, intercept only", {
  
  X <- matrix(0.0, nrow = 100L, ncol = 0L)
  par <- 0.25
  expected <- exp(0.25)*100.0 - 100
  expect_equal(.calEqsGMM(par, X, X), drop(crossprod(expected)))
  
})

test_that("`.calEqGMM()` returns expected results, 1 covariate", {
  
  X <- withr::with_seed(1234L, matrix(stats::rnorm(100), ncol = 1L))
  par <- c(0.25, 0.5)
  wgt <- drop(exp(0.25 + X * 0.5))
  expected <- c(sum(wgt), sum(wgt * X)) - c(100, sum(X))
  expect_equal(.calEqsGMM(par, X, X), drop(crossprod(expected)))
  
})

test_that("`.calEqGMM()` returns expected results, >1 covariate", {
  
  X <- withr::with_seed(1234L, matrix(stats::rnorm(100), ncol = 2L))
  par <- c(0.25, 0.5, 0.75)
  wgt <- drop(exp(0.25 + 0.5 * X[, 1L] + 0.75 * X[, 2L]))
  expected <- c(sum(wgt), sum(wgt * X[, 1L]), sum(wgt * X[, 2L])) - 
    c(50, sum(X[, 1L]), sum(X[, 2L]))
  expect_equal(.calEqsGMM(par, X, X), drop(crossprod(expected)))
  
})

test_that("`.calEqGMM()` returns expected results, >1 covariate", {
  
  X <- withr::with_seed(1234L, matrix(stats::rnorm(100), ncol = 2L))
  par <- c(0.25, 0.5, 0.75)
  wgt <- drop(exp(0.25 + 0.5 * X[, 1L] + 0.75 * X[, 2L]))
  expected <- c(sum(wgt), sum(wgt * X[, 1L]), sum(wgt * X[, 2L])) - 
    c(50, 2.0*sum(X[, 1L]), 2.0*sum(X[, 2L]))
  expect_equal(.calEqsGMM(par, 2.0*X, X), drop(crossprod(expected)))
  
})