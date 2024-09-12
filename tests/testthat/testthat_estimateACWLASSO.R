test_that("`.estimateACWLASSO()` returns expected results multiple EC", {
  
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(500), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "mainName" = paste0("X", 1:3),
      "Y.hat.A0" = rnorm(100, 1.4, 8),
      "Y.hat.A1" = rnorm(100, 1.3, 7)
    )
  )
  data.ec <- list()
  data.ec[[1L]] <- withr::with_seed(
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
  
  data.ec[[2L]] <- withr::with_seed(
    4563, 
    list(
      "X" = matrix(rnorm(10000), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(2000, 1, 8),
      "ps" = runif(2000),
      "mainName" = paste0("X", 1:3),
      "Y.hat" = list("rct" = rnorm(2000, 1.4, 8),
                     "ec" = rnorm(2000, 1.3, 7))
    )
  )
  
  bias <- withr::with_seed(
    2340,
    list(c(rep(0, 500), rep(0.5, 500)),
         c(rep(0, 500), rep(0.5, 1500)))
  )
  
  aipw_result <- .estimateAIPW(data.rct)
  
  acw_result <- .estimateACW(data.rct = data.rct,
                             data.ec = data.ec,
                             aipw.result = aipw_result,
                             bias = bias)
  n_rct <- 100
  
  tau_hat <- numeric(2L)
  tau_hat[1L] <- sum(acw_result$tau.i[[1L]] * c(rep(1, 100), rep(1, 500), rep(0, 500))) / 100
  tau_hat[2L] <- sum(acw_result$tau.i[[2L]] * c(rep(1, 100), rep(1, 500), rep(0, 1500))) / 100
  
  expected <- list("tau.hat" = tau_hat,
                   "sd.hat" = abs(acw_result$sd.hat),
                   "tau.score" = acw_result$tau.score)
  
  expect_equal(.estimateACWLASSO(data.rct, data.ec, aipw_result, bias), expected)
  
})

test_that("`.estimateACWLASSO()` returns expected results 1 EC", {
  
  data.rct <- withr::with_seed(
    1234, 
    list(
      "X" = matrix(rnorm(500), ncol = 5, dimnames = list(NULL, paste0("X", 1:5))),
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "mainName" = paste0("X", 1:3),
      "Y.hat.A0" = rnorm(100, 1.4, 8),
      "Y.hat.A1" = rnorm(100, 1.3, 7)
    )
  )
  data.ec <- list()
  data.ec[[1L]] <- withr::with_seed(
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
  
  bias <- withr::with_seed(
    2340,
    list(c(rep(0, 500), rep(0.5, 500)))
  )
  
  aipw_result <- .estimateAIPW(data.rct)
  
  acw_result <- .estimateACW(data.rct = data.rct,
                             data.ec = data.ec,
                             aipw.result = aipw_result,
                             bias = bias)
  n_rct <- 100
  
  tau_hat <- numeric(1L)
  tau_hat[1L] <- sum(acw_result$tau.i[[1L]] * c(rep(1, 100), rep(1, 500), rep(0, 500))) / 100

  expected <- list("tau.hat" = tau_hat,
                   "sd.hat" = abs(acw_result$sd.hat),
                   "tau.score" = acw_result$tau.score)
  
  expect_equal(.estimateACWLASSO(data.rct, data.ec, aipw_result, bias), expected)
  
})