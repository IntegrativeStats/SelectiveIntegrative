test_that("`estimateACW()` returns expected results multiple EC no bias", {
  
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
  
  aipw_result <- .estimateAIPW(data.rct)
  
  ec_result <- list()
  ec_result[[1L]] <- .estimateEC(data.rct, data.ec[[1L]], NULL)
  ec_result[[2L]] <- .estimateEC(data.rct, data.ec[[2L]], NULL)

  n_ei <- c(1000, 2000)  

  sd_hat <- numeric(2)
  tau_score <- list()
  tau <- list()
  tau_hat <- numeric(2)
  
  mu <- c(aipw_result$mu$A1, rep(0.0, 1000))
  tau[[1L]] <- mu - ec_result[[1L]]$mu0.i
  tau_score[[1L]] <- mu - ec_result[[1L]]$mu0.score
  tau_hat[1L] <- sum(tau[[1L]], na.rm = TRUE) / 100
  sd_hat[1L] <- {sum({tau_score[[1L]] - 
      c(rep(tau_hat[1L], 100), rep(0.0, 1000))}^2, 
      na.rm = TRUE) / 100} |> sqrt()
  
  mu <- c(aipw_result$mu$A1, rep(0.0, 2000))
  tau[[2L]] <- mu - ec_result[[2L]]$mu0.i
  tau_score[[2L]] <- mu - ec_result[[2L]]$mu0.score
  tau_hat[2L] <- sum(tau[[2L]], na.rm = TRUE) / 100
  sd_hat[2L] <- {sum({tau_score[[2L]] - 
      c(rep(tau_hat[2L], 100), rep(0.0, 2000))}^2, 
      na.rm = TRUE) / 100} |> sqrt()
  
  expected <- list("est" = ec_result, 
                   "tau.i" = tau,
                   "tau.hat" = tau_hat,
                   "tau.score" = tau_score,
                   "sd.hat" = sd_hat)
  
  expect_equal(.estimateACW(data.rct, data.ec, aipw_result, NULL), expected)
})

test_that("`estimateACW()` returns expected results multiple EC w/ bias", {
  
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
  
  ec_result <- list()
  ec_result[[1L]] <- .estimateEC(data.rct, data.ec[[1L]], bias[[1L]])
  ec_result[[2L]] <- .estimateEC(data.rct, data.ec[[2L]], bias[[2L]])
  
  n_ei <- c(1000, 2000)  
  
  sd_hat <- numeric(2)
  tau_score <- list()
  tau <- list()
  tau_hat <- numeric(2)
  
  mu <- c(aipw_result$mu$A1, rep(0.0, 1000))
  tau[[1L]] <- mu - ec_result[[1L]]$mu0.i
  tau_score[[1L]] <- mu - ec_result[[1L]]$mu0.score
  tau_hat[1L] <- sum(tau[[1L]], na.rm = TRUE) / 100
  sd_hat[1L] <- {sum({tau_score[[1L]] - 
      c(rep(tau_hat[1L], 100), rep(0.0, 1000))}^2, 
      na.rm = TRUE) / 100} |> sqrt()
  
  mu <- c(aipw_result$mu$A1, rep(0.0, 2000))
  tau[[2L]] <- mu - ec_result[[2L]]$mu0.i
  tau_score[[2L]] <- mu - ec_result[[2L]]$mu0.score
  tau_hat[2L] <- sum(tau[[2L]], na.rm = TRUE) / 100
  sd_hat[2L] <- {sum({tau_score[[2L]] - 
      c(rep(tau_hat[2L], 100), rep(0.0, 2000))}^2, 
      na.rm = TRUE) / 100} |> sqrt()
  
  expected <- list("est" = ec_result, 
                   "tau.i" = tau,
                   "tau.hat" = tau_hat,
                   "tau.score" = tau_score,
                   "sd.hat" = sd_hat)
  
  expect_equal(.estimateACW(data.rct, data.ec, aipw_result, bias), expected)
})

test_that("`estimateACW()` returns expected results 1 EC no bias", {
  
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
  
  aipw_result <- .estimateAIPW(data.rct)
  
  ec_result <- list()
  ec_result[[1L]] <- .estimateEC(data.rct, data.ec[[1L]], NULL)

  n_ei <- c(1000)  
  
  sd_hat <- numeric(1)
  tau_score <- list()
  tau <- list()
  tau_hat <- numeric(1)
  
  mu <- c(aipw_result$mu$A1, rep(0.0, 1000))
  tau[[1L]] <- mu - ec_result[[1L]]$mu0.i
  tau_score[[1L]] <- mu - ec_result[[1L]]$mu0.score
  tau_hat[1L] <- sum(tau[[1L]], na.rm = TRUE) / 100
  sd_hat[1L] <- {sum({tau_score[[1L]] - 
      c(rep(tau_hat[1L], 100), rep(0.0, 1000))}^2, 
      na.rm = TRUE) / 100} |> sqrt()
  
  expected <- list("est" = ec_result, 
                   "tau.i" = tau,
                   "tau.hat" = tau_hat,
                   "tau.score" = tau_score,
                   "sd.hat" = sd_hat)
  
  expect_equal(.estimateACW(data.rct, data.ec, aipw_result, NULL), expected)
})

test_that("`estimateACW()` returns expected results 1 ec w/ bias", {
  
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
  
  ec_result <- list()
  ec_result[[1L]] <- .estimateEC(data.rct, data.ec[[1L]], bias[[1L]])

  n_ei <- c(1000)  
  
  sd_hat <- numeric(1)
  tau_score <- list()
  tau <- list()
  tau_hat <- numeric(1)
  
  mu <- c(aipw_result$mu$A1, rep(0.0, 1000))
  tau[[1L]] <- mu - ec_result[[1L]]$mu0.i
  tau_score[[1L]] <- mu - ec_result[[1L]]$mu0.score
  tau_hat[1L] <- sum(tau[[1L]], na.rm = TRUE) / 100
  sd_hat[1L] <- {sum({tau_score[[1L]] - 
      c(rep(tau_hat[1L], 100), rep(0.0, 1000))}^2, 
      na.rm = TRUE) / 100} |> sqrt()
  
  expected <- list("est" = ec_result, 
                   "tau.i" = tau,
                   "tau.hat" = tau_hat,
                   "tau.score" = tau_score,
                   "sd.hat" = sd_hat)
  
  expect_equal(.estimateACW(data.rct, data.ec, aipw_result, bias), expected)
})