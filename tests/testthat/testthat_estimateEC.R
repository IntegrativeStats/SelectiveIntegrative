test_that("`.lambdahat()` returns expected results", {
  
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(5000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "mainName" = c("X1", "X2", "X3", "X4")
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(rnorm(10000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "mainName" = c("X1", "X2", "X3", "X4")
    )
  )

  p_main_effect <- 5L  
  n_ec <- 10000

  lambda_hat <- nleqslv::nleqslv(x = rep(0.0, 5),
                                 fn = .calEqs,
                                 jac = .calEqsGradient,
                                 X.rct = data.rct$X[, c("X1", "X2", "X3", "X4")], 
                                 X.ec = data.ec$X[, c("X1", "X2", "X3", "X4")], 
                                 control = list(trace = FALSE,
                                                maxit = 500L),
                                 method = 'Newton')$x
  
  expect_equal(.lambdaHat(data.rct, data.ec), lambda_hat)
})

test_that("`.lambdahat()` returns expected results intercept only models", {
  
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(5000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "mainName" = NULL
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(rnorm(10000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "mainName" = NULL
    )
  )
  
  p_main_effect <- 1 
  n_ec <- 10000
  
  lambda_hat <- BB::BBoptim(par = 0.0,
                            fn = .calEqs,
                            X.rct = matrix(0, nrow(data.rct$X), ncol = 0), 
                            X.ec = matrix(0, nrow(data.ec$X), ncol = 0), 
                            control = list(trace = TRUE,
                                           ftol = 1e-8,
                                           maxit = 500L))$par

  expect_equal(.lambdaHat(data.rct, data.ec), lambda_hat)
})

test_that("`.lambdahat()` returns expected results no covariates", {
  
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(0, nrow = 1000L, ncol = 0L),
      "mainName" = NULL
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(0, nrow = 10000L, ncol = 0L),
      "mainName" = NULL
    )
  )
  
  p_main_effect <- 1 
  n_ec <- 10000
  
  lambda_hat <- BB::BBoptim(par = 0.0,
                            fn = .calEqs,
                            X.rct = matrix(0, 1000, ncol = 0), 
                            X.ec = matrix(0, 10000, ncol = 0), 
                            control = list(trace = TRUE,
                                           ftol = 1e-8,
                                           maxit = 500L))$par
  
  expect_equal(.lambdaHat(data.rct, data.ec), lambda_hat)
})

test_that("`.lambdahat()` returns expected results", {
  
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(500), ncol = 50, dimnames = list(NULL, paste0("X", 1:50))),
      "mainName" = paste0("X", 1:50)
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(rnorm(500), ncol = 50, dimnames = list(NULL, paste0("X", 1:50))),
      "mainName" = paste0("X", 1:50)
    )
  )
  
  p_main_effect <- 51L
  n_ec <- 10
  
  lambda_hat <- stats::optim(par = rep(0.0, 51L),
                             fn = .calEqsGMM,
                             X.rct = data.rct$X, 
                             X.ec = data.ec$X, 
                             control = list(trace = TRUE,
                                            abstol = 1e-8,
                                            maxit = 500L),
                             method = 'BFGS')$par

  expect_equal(.lambdaHat(data.rct, data.ec), lambda_hat)
})

.SACW <- function(X.rct, X.ec, q.hat.ec) {
  # {n_rct + n_ec x p_me}
  S_acw_q <- rbind(-X.rct, X.ec * q.hat.ec)
  # {p_me x p_me}
  dot_S_acw_q <- t(q.hat.ec * X.ec) %*% X.ec / nrow(X.rct)
  
  # {p_me x p_me}
  inv_inf <- tryCatch(MASS::ginv(dot_S_acw_q),
                      error = function(e) {
                        stop("error encountered in evaluating inverse of influence function\n\t",
                             e$message, call. = FALSE)
                      })
  
  list("S" = S_acw_q, "inv.inf" = inv_inf)
  
}

test_that("`.SACW()` returns expected results", {
  
  X.rct <- withr::with_seed(
    1234, matrix(rnorm(500), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))))
  
  X.ec <- withr::with_seed(
    2456, matrix(rnorm(5000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))))
  
  q.hat.ec <- withr::with_seed(456, rnorm(1000))
  
  S_acw_q <- matrix(0.0, 1100, 5, dimnames = list(NULL, paste0("X", 1:5)))
  S_acw_q[1:100, 1:5] <- -X.rct
  S_acw_q[101:1100, 1:5] <- cbind(X.ec[,1]*q.hat.ec,
                                  X.ec[,2]*q.hat.ec,
                                  X.ec[,3]*q.hat.ec,
                                  X.ec[,4]*q.hat.ec,
                                  X.ec[,5]*q.hat.ec)
  dot_S_acw_q <- rbind(X.ec[,1]*q.hat.ec,
                       X.ec[,2]*q.hat.ec,
                       X.ec[,3]*q.hat.ec,
                       X.ec[,4]*q.hat.ec,
                       X.ec[,5]*q.hat.ec) %*% X.ec / 100.0
  
  expected <- list("S" = S_acw_q, "inv.inf" = MASS::ginv(dot_S_acw_q))

  expect_equal(.SACW(X.rct, X.ec, q.hat.ec), expected)
  
})

test_that("`.SACW()` returns expected results intercept only", {
  
  X.rct <- matrix(1.0, nrow = 100L, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  
  X.ec <- matrix(1.0, nrow = 1000L, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  
  q.hat.ec <- withr::with_seed(456, rnorm(1000))
  
  S_acw_q <- matrix(0.0, 1100, 1, dimnames = list(NULL, "(Intercept)"))
  S_acw_q[1:100, 1L] <- -1.0
  S_acw_q[101:1100, 1L] <- q.hat.ec
  
  dot_S_acw_q <- {X.ec[,1]*q.hat.ec} %*% X.ec / 100.0
  
  expected <- list("S" = S_acw_q, "inv.inf" = MASS::ginv(dot_S_acw_q))

  expect_equal(.SACW(X.rct, X.ec, q.hat.ec), expected)
  
})

test_that("`.influenceFunction()` returns expected errors", {
  
  expect_error(.influenceFunction(),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  expect_error(.influenceFunction(data.frame("X" = 1, "A" = 1, "Y" = 1)),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct <- list()
  expect_error(.influenceFunction(data.rct),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct <- list("X" = matrix(1, 10, 3, dimnames = list(NULL, c("x1", "x2", "x3"))))
  expect_error(.influenceFunction(data.rct),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$A <- rep(1, 10)
  expect_error(.influenceFunction(data.rct),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$Y <- rep(1, 10)
  expect_error(.influenceFunction(data.rct),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$mainName <- c("x1", "x2")
  expect_error(.influenceFunction(data.rct),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$ps <- rep(1, 10)
  expect_error(.influenceFunction(data.rct),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$Y.hat <- rep(1, 10)
  expect_error(.influenceFunction(data.rct),
               "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$Y.hat <- NULL
  data.rct$Y.hat.A0 <- rep(1, 10)
  

  expect_error(.influenceFunction(data.rct),
               "`data.ec` must be a list containing {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  expect_error(.influenceFunction(data.rct, data.frame("X" = 1, "A" = 1, "Y" = 1)),
               "`data.ec` must be a list containing {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec <- list()
  expect_error(.influenceFunction(data.rct, data.ec),
               "`data.ec` must be a list containing {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec <- list("X" = matrix(1, 10, 3, dimnames = list(NULL, c("x1", "x2", "x3"))))
  expect_error(.influenceFunction(data.rct, data.ec),
               "`data.ec` must be a list containing {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec$Y <- rep(1, 10)
  expect_error(.influenceFunction(data.rct, data.ec),
               "`data.ec` must be a list containing {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec$mainName <- c("x1", "x2")
  expect_error(.influenceFunction(data.rct, data.ec),
               "`data.ec` must be a list containing {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec$ps <- rep(1, 10)
  expect_error(.influenceFunction(data.rct, data.ec),
               "`data.ec` must be a list containing {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec$Y.hat <- rep(1, 10)

  expect_error(.influenceFunction(data.rct, data.ec),
               "`data.ec$Y.hat` must be a list containing element {rct}",
               fixed = TRUE)
  data.ec$Y.hat <- list(1:10)
  expect_error(.influenceFunction(data.rct, data.ec),
               "`data.ec$Y.hat` must be a list containing element {rct}",
               fixed = TRUE)
  data.ec$Y.hat <- list("RCT" = 1:10)
  expect_error(.influenceFunction(data.rct, data.ec),
               "`data.ec$Y.hat` must be a list containing element {rct}",
               fixed = TRUE)
  data.ec$Y.hat <- list("rct" = rep(1, 10))
  
  expect_error(.influenceFunction(data.rct, data.ec), 
               "`lambda.hat` must be a numeric vector of lenght p_main_effects")
  lambda.hat <- c("a", "b", "c")
  expect_error(.influenceFunction(data.rct, data.ec, lambda.hat), 
               "`lambda.hat` must be a numeric vector of lenght p_main_effects")
  lambda.hat <- 1:2
  expect_error(.influenceFunction(data.rct, data.ec, lambda.hat), 
               "`lambda.hat` must be a numeric vector of lenght p_main_effects")
  lambda.hat <- 1:3

  expect_error(.influenceFunction(data.rct, data.ec, lambda.hat), 
               "`r.X` must be a positive scalar numeric")
  r.X <- 1:2
  expect_error(.influenceFunction(data.rct, data.ec, lambda.hat), 
               "`r.X` must be a positive scalar numeric")
  r.X <- -0.1
  expect_error(.influenceFunction(data.rct, data.ec, lambda.hat), 
               "`r.X` must be a positive scalar numeric")
  r.X <- 1e3
  
  expect_error(.influenceFunction(data.rct, data.ec, lambda.hat, r.X), 
               "`zero.bias` must be a logical vector")
  zero.bias <- rep(1, 10)
  expect_error(.influenceFunction(data.rct, data.ec, lambda.hat, r.X, zero.bias), 
               "`zero.bias` must be a logical vector")
  zero.bias <- rep(TRUE, 9)
  expect_error(.influenceFunction(data.rct, data.ec, lambda.hat, r.X, zero.bias), 
               "`zero.bias` must be a logical vector")
})

test_that(".influenceFunction() returns expected results multiple covariates", {
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(500), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "mainName" = paste0("X", 1:3),
      "Y.hat.A0" = rnorm(100, 1.4, 8)
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(rnorm(5000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(1000, 1, 8),
      "ps" = runif(1000),
      "mainName" = paste0("X", 1:3),
      "Y.hat" = list("rct" = rnorm(1000, 1.4, 8))
    )
  )

  lambda.hat <- c(0.25, 0.6, -0.8, 1.25)  
  r.X <- 0.1
  zero.bias <- withr::with_seed(321,
                                as.logical(rbinom(1000, 1, 0.8)))
  
  q_hat_ec <- exp(lambda.hat[1L] + lambda.hat[2L] * data.ec$X[, "X1"] +
                  lambda.hat[3L] * data.ec$X[, "X2"] +
                  lambda.hat[4L] * data.ec$X[, "X3"])
  q_hat_ec[q_hat_ec > 1e8] <- 1e8
  
  q_hat_rct <- exp(lambda.hat[1L] + lambda.hat[2L] * data.rct$X[, "X1"] +
                     lambda.hat[3L] * data.rct$X[, "X2"] +
                     lambda.hat[4L] * data.rct$X[, "X3"])
  q_hat_rct[q_hat_rct > 1e8] <- 1e8
  
  ME_rct_with_intercept <- cbind(1.0, data.rct$X[, "X1"], data.rct$X[, "X2"], data.rct$X[, "X3"])
  ME_ec_with_intercept <- cbind(1.0, data.ec$X[, "X1"], data.ec$X[, "X2"], data.ec$X[, "X3"])
  
  S_results <- .SACW(X.rct = ME_rct_with_intercept, 
                     X.ec = ME_ec_with_intercept, 
                     q.hat.ec = q_hat_ec)
    
  deno_rct <- r.X + {1.0 - data.rct$ps} * q_hat_rct
  numo_rct <- q_hat_rct *  {1.0 - data.rct$A} * {data.rct$Y - data.rct$Y.hat.A0}
  tmp <- r.X * numo_rct / deno_rct^2 * ME_rct_with_intercept
  dot_mu0_rct_q <- colSums(tmp, na.rm = TRUE) / 100.0
    
  deno_ec <- r.X * zero.bias + {1.0 - data.ec$ps} * q_hat_ec
  numo_ec <- q_hat_ec * r.X * zero.bias * {data.ec$Y - data.ec$Y.hat$rct}
  tmp <- r.X * numo_ec / deno_ec^2 * ME_ec_with_intercept
  dot_mu0_ec_q <- colSums(tmp, na.rm = TRUE) / 100.0
    
  mu0_i <- c(data.rct$Y.hat.A0 + numo_rct / deno_rct, numo_ec / deno_ec)
    
  mu0_score <-  mu0_i - drop(S_results$S %*% t(dot_mu0_rct_q %*% S_results$inv.inf)) -
    drop(S_results$S %*% t(dot_mu0_ec_q %*% S_results$inv.inf))
   
  expected <-   list("q.hat.ec" = q_hat_ec,
                     "mu0.i" = mu0_i,
                     "mu0.score" = mu0_score)
    
  expect_equal(.influenceFunction(data.rct, data.ec, lambda.hat, r.X, zero.bias),
               expected)
  
})

test_that(".influenceFunction() returns expected results single covariates", {
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(100), ncol = 1, dimnames = list(NULL, paste0("X", 1))),
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "mainName" = "X1",
      "Y.hat.A0" = rnorm(100, 1.4, 8)
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(rnorm(1000), ncol = 1, dimnames = list(NULL, paste0("X", 1))),
      "Y" = rnorm(1000, 1, 8),
      "ps" = runif(1000),
      "mainName" = "X1",
      "Y.hat" = list("rct" = rnorm(1000, 1.4, 8))
    )
  )
  
  lambda.hat <- c(0.25, 0.6)  
  r.X <- 0.1
  zero.bias <- withr::with_seed(321,
                                as.logical(rbinom(1000, 1, 0.8)))
  
  q_hat_ec <- exp(lambda.hat[1L] + lambda.hat[2L] * data.ec$X[, "X1"])
  q_hat_ec[q_hat_ec > 1e8] <- 1e8
  
  q_hat_rct <- exp(lambda.hat[1L] + lambda.hat[2L] * data.rct$X[, "X1"])
  q_hat_rct[q_hat_rct > 1e8] <- 1e8
  
  ME_rct_with_intercept <- cbind(1.0, data.rct$X[, "X1"])
  ME_ec_with_intercept <- cbind(1.0, data.ec$X[, "X1"])
  
  S_results <- .SACW(X.rct = ME_rct_with_intercept, 
                     X.ec = ME_ec_with_intercept, 
                     q.hat.ec = q_hat_ec)
  
  deno_rct <- r.X + {1.0 - data.rct$ps} * q_hat_rct
  numo_rct <- q_hat_rct *  {1.0 - data.rct$A} * {data.rct$Y - data.rct$Y.hat.A0}
  tmp <- r.X * numo_rct / deno_rct^2 * ME_rct_with_intercept
  dot_mu0_rct_q <- colSums(tmp, na.rm = TRUE) / 100.0
  
  deno_ec <- r.X * zero.bias + {1.0 - data.ec$ps} * q_hat_ec
  numo_ec <- q_hat_ec * r.X * zero.bias * {data.ec$Y - data.ec$Y.hat$rct}
  tmp <- r.X * numo_ec / deno_ec^2 * ME_ec_with_intercept
  dot_mu0_ec_q <- colSums(tmp, na.rm = TRUE) / 100.0
  
  mu0_i <- c(data.rct$Y.hat.A0 + numo_rct / deno_rct, numo_ec / deno_ec)
  
  mu0_score <-  mu0_i - drop(S_results$S %*% t(dot_mu0_rct_q %*% S_results$inv.inf)) -
    drop(S_results$S %*% t(dot_mu0_ec_q %*% S_results$inv.inf))
  
  expected <-   list("q.hat.ec" = q_hat_ec,
                     "mu0.i" = mu0_i,
                     "mu0.score" = mu0_score)
  
  expect_equal(.influenceFunction(data.rct, data.ec, lambda.hat, r.X, zero.bias),
               expected)
  
})

test_that(".influenceFunction() returns expected results intercept model single covariates", {
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(100), ncol = 1, dimnames = list(NULL, paste0("X", 1))),
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "mainName" = NULL,
      "Y.hat.A0" = rnorm(100, 1.4, 8)
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(rnorm(1000), ncol = 1, dimnames = list(NULL, paste0("X", 1))),
      "Y" = rnorm(1000, 1, 8),
      "ps" = runif(1000),
      "mainName" = NULL,
      "Y.hat" = list("rct" = rnorm(1000, 1.4, 8))
    )
  )
  
  lambda.hat <- c(0.25)  
  r.X <- 0.1
  zero.bias <- withr::with_seed(321,
                                as.logical(rbinom(1000, 1, 0.8)))
  
  q_hat_ec <- exp(lambda.hat[1L]) |> rep(1000)
  q_hat_ec[q_hat_ec > 1e8] <- 1e8
  
  q_hat_rct <- exp(lambda.hat[1L]) |> rep(100)
  q_hat_rct[q_hat_rct > 1e8] <- 1e8
  
  ME_rct_with_intercept <- matrix(1.0, nrow = 100, ncol = 1L)
  ME_ec_with_intercept <- matrix(1.0, nrow = 1000, ncol = 1L)
  
  S_results <- .SACW(X.rct = ME_rct_with_intercept, 
                     X.ec = ME_ec_with_intercept, 
                     q.hat.ec = q_hat_ec)
  
  deno_rct <- r.X + {1.0 - data.rct$ps} * q_hat_rct
  numo_rct <- q_hat_rct *  {1.0 - data.rct$A} * {data.rct$Y - data.rct$Y.hat.A0}
  tmp <- r.X * numo_rct / deno_rct^2 * ME_rct_with_intercept
  dot_mu0_rct_q <- colSums(tmp, na.rm = TRUE) / 100.0
  
  deno_ec <- r.X * zero.bias + {1.0 - data.ec$ps} * q_hat_ec
  numo_ec <- q_hat_ec * r.X * zero.bias * {data.ec$Y - data.ec$Y.hat$rct}
  tmp <- r.X * numo_ec / deno_ec^2 * ME_ec_with_intercept
  dot_mu0_ec_q <- colSums(tmp, na.rm = TRUE) / 100.0
  
  mu0_i <- c(data.rct$Y.hat.A0 + numo_rct / deno_rct, numo_ec / deno_ec)
  
  mu0_score <-  mu0_i - drop(S_results$S %*% t(dot_mu0_rct_q %*% S_results$inv.inf)) -
    drop(S_results$S %*% t(dot_mu0_ec_q %*% S_results$inv.inf))
  
  expected <-   list("q.hat.ec" = q_hat_ec,
                     "mu0.i" = mu0_i,
                     "mu0.score" = mu0_score)
  
  expect_equal(.influenceFunction(data.rct, data.ec, lambda.hat, r.X, zero.bias),
               expected)
  
})

test_that(".influenceFunction() returns expected results intercept model no covariates", {
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(0, nrow = 100, ncol = 0),
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "mainName" = NULL,
      "Y.hat.A0" = rnorm(100, 1.4, 8)
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(0, nrow = 1000, ncol = 0),
      "Y" = rnorm(1000, 1, 8),
      "ps" = runif(1000),
      "mainName" = NULL,
      "Y.hat" = list("rct" = rnorm(1000, 1.4, 8))
    )
  )
  
  lambda.hat <- c(0.25)  
  r.X <- 0.1
  zero.bias <- withr::with_seed(321,
                                as.logical(rbinom(1000, 1, 0.8)))
  
  q_hat_ec <- exp(lambda.hat[1L]) |> rep(1000)
  q_hat_ec[q_hat_ec > 1e8] <- 1e8
  
  q_hat_rct <- exp(lambda.hat[1L]) |> rep(100)
  q_hat_rct[q_hat_rct > 1e8] <- 1e8
  
  ME_rct_with_intercept <- matrix(1.0, nrow = 100, ncol = 1L)
  ME_ec_with_intercept <- matrix(1.0, nrow = 1000, ncol = 1L)
  
  S_results <- .SACW(X.rct = ME_rct_with_intercept, 
                     X.ec = ME_ec_with_intercept, 
                     q.hat.ec = q_hat_ec)
  
  deno_rct <- r.X + {1.0 - data.rct$ps} * q_hat_rct
  numo_rct <- q_hat_rct *  {1.0 - data.rct$A} * {data.rct$Y - data.rct$Y.hat.A0}
  tmp <- r.X * numo_rct / deno_rct^2 * ME_rct_with_intercept
  dot_mu0_rct_q <- colSums(tmp, na.rm = TRUE) / 100.0
  
  deno_ec <- r.X * zero.bias + {1.0 - data.ec$ps} * q_hat_ec
  numo_ec <- q_hat_ec * r.X * zero.bias * {data.ec$Y - data.ec$Y.hat$rct}
  tmp <- r.X * numo_ec / deno_ec^2 * ME_ec_with_intercept
  dot_mu0_ec_q <- colSums(tmp, na.rm = TRUE) / 100.0
  
  mu0_i <- c(data.rct$Y.hat.A0 + numo_rct / deno_rct, numo_ec / deno_ec)
  
  mu0_score <-  mu0_i - drop(S_results$S %*% t(dot_mu0_rct_q %*% S_results$inv.inf)) -
    drop(S_results$S %*% t(dot_mu0_ec_q %*% S_results$inv.inf))
  
  expected <-   list("q.hat.ec" = q_hat_ec,
                     "mu0.i" = mu0_i,
                     "mu0.score" = mu0_score)
  
  expect_equal(.influenceFunction(data.rct, data.ec, lambda.hat, r.X, zero.bias),
               expected)
  
})

test_that("`.estimateEC()` returns expected errors", {
  
  expect_error(.estimateEC(),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  expect_error(.estimateEC(data.frame("X" = 1, "A" = 1, "Y" = 1)),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct <- list()
  expect_error(.estimateEC(data.rct),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct <- list("X" = matrix(1, 10, 3, dimnames = list(NULL, c("x1", "x2", "x3"))))
  expect_error(.estimateEC(data.rct),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$A <- rep(1, 10)
  expect_error(.estimateEC(data.rct),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$Y <- rep(1, 10)
  expect_error(.estimateEC(data.rct),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$mainName <- c("x1", "x2")
  expect_error(.estimateEC(data.rct),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$ps <- rep(1, 10)
  expect_error(.estimateEC(data.rct),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$Y.hat <- rep(1, 10)
  expect_error(.estimateEC(data.rct),
               "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}",
               fixed = TRUE)
  data.rct$Y.hat <- NULL
  data.rct$Y.hat.A0 <- rep(1, 10)
  
  
  expect_error(.estimateEC(data.rct),
               "`data.ec` must be a list containing elements {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  expect_error(.estimateEC(data.rct, data.frame("X" = 1, "A" = 1, "Y" = 1)),
               "`data.ec` must be a list containing elements {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec <- list()
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec` must be a list containing elements {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec <- list("X" = matrix(1, 10, 3, dimnames = list(NULL, c("x1", "x2", "x3"))))
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec` must be a list containing elements {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec$Y <- rep(1, 10)
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec` must be a list containing elements {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec$mainName <- c("x1", "x2")
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec` must be a list containing elements {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec$ps <- rep(1, 10)
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec` must be a list containing elements {X, Y, mainName, ps, Y.hat}",
               fixed = TRUE)
  data.ec$Y.hat <- rep(1, 10)
  
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec$Y.hat` must be a list containing elements {rct, ec}",
               fixed = TRUE)
  data.ec$Y.hat <- list(1:10)
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec$Y.hat` must be a list containing elements {rct, ec}",
               fixed = TRUE)
  data.ec$Y.hat <- list("RCT" = 1:10)
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec$Y.hat` must be a list containing elements {rct, ec}",
               fixed = TRUE)
  data.ec$Y.hat <- list("rct" = rep(1, 10))
  expect_error(.estimateEC(data.rct, data.ec),
               "`data.ec$Y.hat` must be a list containing elements {rct, ec}",
               fixed = TRUE)
  data.ec$Y.hat$ec <- 1:10
  

  expect_error(.estimateEC(data.rct, data.ec, NA),
               "`bias` must be NULL or a numeric vector")  
  expect_error(.estimateEC(data.rct, data.ec, rep(TRUE, 10)),
               "`bias` must be NULL or a numeric vector") 
  expect_error(.estimateEC(data.rct, data.ec, rep(1, 9)),
               "`bias` must be NULL or a numeric vector") 
 
})

test_that("`.estimateEC()` returns expected results multiple covariates", {
  
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(500), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "mainName" = paste0("X", 1:3),
      "Y.hat.A0" = rnorm(100, 1.4, 8)
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(rnorm(5000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(1000, 1, 8),
      "ps" = runif(1000),
      "mainName" = paste0("X", 1:3),
      "Y.hat" = list("rct" = rnorm(1000, 1.4, 8),
                     "ec" = rnorm(1000, 1.3, 7))
    )
  )
  
  lambda_hat <- .lambdaHat(data.rct = data.rct, data.ec = data.ec)
    
  bias <- rep(0.0, 1000)
  zero_bias <- rep(TRUE, 1000)

  tmp <- mean({data.ec$Y.hat$ec - mean(data.ec$Y)}^2)
  r_X <- mean({data.rct$Y.hat.A0[data.rct$A == 0L] - mean(data.rct$Y[data.rct$A == 0L])}^2) / tmp
  
  bias <- data.ec$Y - data.ec$Y.hat$rct
  bias_i_score <- c(rep(0.0, 100), bias)
    
  mu0 <- .influenceFunction(data.rct = data.rct, 
                            data.ec = data.ec, 
                            lambda.hat = lambda_hat, 
                            r.X = r_X, 
                            zero.bias = zero_bias)
    
  expected <-  list("q.hat.ec" = mu0$q.hat.ec,
                    "mu0.i" = mu0$mu0.i,
                    "mu0.score" = mu0$mu0.score,
                    "bias" = bias,
                    "bias.score" = bias_i_score,
                    "gamma.hat" = data.ec$Y.hat$ec - data.ec$Y.hat$rct, 
                    "r.X" = r_X)

  expect_equal(.estimateEC(data.rct, data.ec, NULL), expected)
})

test_that("`.estimateEC()` returns expected results multiple covariates", {
  
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(500), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "mainName" = paste0("X", 1:3),
      "Y.hat.A0" = rnorm(100, 1.4, 8)
    )
  )
  
  data.ec <- withr::with_seed(
    2456, 
    list(
      "X" = matrix(rnorm(5000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(1000, 1, 8),
      "ps" = runif(1000),
      "mainName" = paste0("X", 1:3),
      "Y.hat" = list("rct" = rnorm(1000, 1.4, 8),
                     "ec" = rnorm(1000, 1.3, 7))
    )
  )
  
  lambda_hat <- .lambdaHat(data.rct = data.rct, data.ec = data.ec)
  
  bias_in <- c(rep(0.0, 500), rep(0.5, 500))
  zero_bias <- c(rep(TRUE, 500), rep(FALSE, 500))
  
  tmp <- mean({data.ec$Y.hat$ec[1L:500L] - mean(data.ec$Y[1L:500L])}^2)
  r_X <- mean({data.rct$Y.hat.A0[data.rct$A == 0L] - mean(data.rct$Y[data.rct$A == 0L])}^2) / tmp
  
  bias <- data.ec$Y - data.ec$Y.hat$rct
  bias_i_score <- c(rep(0.0, 100), bias)

  mu0 <- .influenceFunction(data.rct = data.rct, 
                            data.ec = data.ec, 
                            lambda.hat = lambda_hat, 
                            r.X = r_X, 
                            zero.bias = zero_bias)
  
  expected <-  list("q.hat.ec" = mu0$q.hat.ec,
                    "mu0.i" = mu0$mu0.i,
                    "mu0.score" = mu0$mu0.score,
                    "bias" = bias,
                    "bias.score" = bias_i_score,
                    "gamma.hat" = data.ec$Y.hat$ec - data.ec$Y.hat$rct, 
                    "r.X" = r_X)

  expect_equal(.estimateEC(data.rct, data.ec, bias_in), expected)
})