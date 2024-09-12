test_that("`.tradeoff()` returns expected results multiple EC", {
  bias <- withr::with_seed(
    2340,
    list(c(rep(0, 500), rep(0.5, 500)),
         c(rep(0, 500), rep(0.5, 1500)))
  )
  
  bias_hat <- c(1, rep(0, 500), rep(0.35, 500), rep(0, 400), rep(0.2, 1600))
  
  bias_score <- withr::with_seed(
    2340,
    list(c(rnorm(1100, 1, 3)),
         c(rnorm(2100, 2, 1)))
  )
  
  zero_x <- c(rep(TRUE, 500), rep(FALSE, 500), rep(TRUE, 400), rep(FALSE, 1600))
  tau_acw_lin <- c(bias[[1L]] * c(rep(1, 500), rep(0, 500)),
                   bias[[2L]] * c(rep(1, 400), rep(0, 1600)))
    
  bias_high <- c(sum(bias[[1L]] * c(rep(1, 500), rep(0, 500))), 
                 sum(bias[[2L]] * c(rep(1, 400), rep(0, 1600)))) / 100
  bias_high <- sum(bias_high^2)
    
  tau_acw_lin_list <- list("1" = (bias_score[[1L]][101:1100] * c(rep(1, 500), rep(0, 500))),
                           "2" = (bias_score[[2L]][101:2100] * c(rep(1, 400), rep(0, 1600))))

  tau_acw_score_recont <- list("1" = c(bias_score[[1L]][1:100], tau_acw_lin_list[[1L]]),
                               "2" = c(bias_score[[2L]][1:100], tau_acw_lin_list[[2L]]))
  
  var_i <- unlist(lapply(tau_acw_score_recont,
                         function(x) { sum({x - mean(x)}^2) / 100 }))
  expected <- mean(var_i) + bias_high
  
  expect_equal(.tradeOff(bias_hat, bias, bias_score, 100), expected)
})

test_that("`.tradeoff()` returns expected results single EC", {
  bias <- withr::with_seed(
    2340,
    list(c(rep(0, 500), rep(0.5, 500)))
  )
  
  bias_hat <- c(1, rep(0, 500), rep(0.35, 500))
  
  bias_score <- withr::with_seed(
    2340,
    list(c(rnorm(1100, 1, 3)))
  )
  
  zero_x <- c(rep(TRUE, 500), rep(FALSE, 500))
  tau_acw_lin <- c(bias[[1L]] * c(rep(1, 500), rep(0, 500)))
  
  bias_high <- c(sum(bias[[1L]] * c(rep(1, 500), rep(0, 500)))) / 100
  bias_high <- sum(bias_high^2)
  
  tau_acw_lin_list <- list("1" = (bias_score[[1L]][101:1100] * c(rep(1, 500), rep(0, 500))))
  
  tau_acw_score_recont <- list("1" = c(bias_score[[1L]][1:100], tau_acw_lin_list[[1L]]))
  
  var_i <- unlist(lapply(tau_acw_score_recont,
                         function(x) { sum({x - mean(x)}^2) / 100 }))
  expected <- mean(var_i) + bias_high
  
  expect_equal(.tradeOff(bias_hat, bias, bias_score, 100), expected)
})

test_that(".cvLasso() returns expected results multiple EC", {
  bias <- withr::with_seed(
    2340,
    list(c(rep(0, 500), rep(0.5, 500)),
         c(rep(0, 500), rep(0.5, 1500)))
  )
  
  bias_hat <- c(1, rep(0, 500), rep(0.35, 500), rep(0, 400), rep(0.2, 1600))
  
  bias_score <- withr::with_seed(
    2340,
    list(c(rnorm(1100, 1, 3)),
         c(rnorm(2100, 2, 1)))
  )
  
  gamma.hat <- withr::with_seed(34523, rnorm(3000))
  nu <- 1.0
  
  Z_mat <- diag(3000)
  min_lambda <- {3000}^-(1.0 + 1.0 / 5.0) / 2.0
  max_lambda <- {3000}^{-0.1}
  
  lambda_vector <- seq(min_lambda, max_lambda, length.out = 100L)

  b_k_w <- abs(gamma.hat)^{-nu}
  b_k_w <- b_k_w / sum(b_k_w) * 3000
  b_k_w[is.na(b_k_w)] <- Inf
  
  fit_lasso <- withr::with_seed(1234,
                                glmnet::glmnet(Z_mat, unlist(bias), 
                                               family = "gaussian",
                                               alpha = 1,
                                               penalty.factor = b_k_w,
                                               intercept = FALSE,
                                               standardize = FALSE,
                                               standardize.response = FALSE,
                                               lambda = lambda_vector))

  lasso_prediction <- predict(fit_lasso, newx = Z_mat, type = "coef")
  
  trade_off <- apply(lasso_prediction, MARGIN = 2L,
                     .tradeOff,
                     bias = bias, 
                     bias.score = bias_score, 
                     n.rct = 100)

  lambda_opt_idx <- which.min(trade_off)
    
  bias_lasso <- lasso_prediction[-1L, lambda_opt_idx] |> drop()
    
  expected <- list("min" = trade_off[lambda_opt_idx],
                   "bias" = bias_lasso)
  
  expect_equal(withr::with_seed(1234,.cvLASSO(nu, gamma.hat, bias, 100, bias_score, NULL)), expected)
  
  
})

test_that(".cvLasso() returns expected results single EC", {
  bias <- withr::with_seed(
    2340,
    list(c(rep(0, 500), rep(0.5, 500)))
  )
  
  bias_hat <- c(1, rep(0, 500), rep(0.35, 500))
  
  bias_score <- withr::with_seed(
    2340,
    list(c(rnorm(1100, 1, 3)))
  )
  
  gamma.hat <- withr::with_seed(34523, rnorm(1000))
  nu <- 1.0
  
  Z_mat <- diag(1000)
  min_lambda <- {1000}^-(1.0 + 1.0 / 5.0) / 2.0
  max_lambda <- {1000}^{-0.1}
  
  lambda_vector <- seq(min_lambda, max_lambda, length.out = 100L)
  
  b_k_w <- abs(gamma.hat)^{-nu}
  b_k_w <- b_k_w / sum(b_k_w) * 1000
  b_k_w[is.na(b_k_w)] <- Inf
  
  fit_lasso <- withr::with_seed(1234,
                                glmnet::glmnet(Z_mat, unlist(bias), 
                                               family = "gaussian",
                                               alpha = 1,
                                               penalty.factor = b_k_w,
                                               intercept = FALSE,
                                               standardize = FALSE,
                                               standardize.response = FALSE,
                                               lambda = lambda_vector))
  
  lasso_prediction <- predict(fit_lasso, newx = Z_mat, type = "coef")
  
  trade_off <- apply(lasso_prediction, MARGIN = 2L,
                     .tradeOff,
                     bias = bias, 
                     bias.score = bias_score, 
                     n.rct = 100)
  
  lambda_opt_idx <- which.min(trade_off)
  
  bias_lasso <- lasso_prediction[-1L, lambda_opt_idx] |> drop()
  
  expected <- list("min" = trade_off[lambda_opt_idx],
                   "bias" = bias_lasso)
  
  expect_equal(withr::with_seed(1234,.cvLASSO(nu, gamma.hat, bias, 100, bias_score, NULL)), expected)
  
  
})