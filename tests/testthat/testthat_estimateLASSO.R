test_that("`.estimateLASSO()` returns expected results multiple EC", {
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
  
  acw_result <- .estimateACW(data.rct, data.ec, aipw_result, bias)
  
  bias_ec <- list(
    acw_result$est[[1L]]$bias,
    acw_result$est[[2L]]$bias
  )

  n_ei <- c(1000, 2000)
  K <- 2L

  bias_ec_score <- list(
    acw_result$est[[1L]]$bias.score,
    acw_result$est[[2L]]$bias.score
  )
    
  gamma_hat <- c(acw_result$est[[1L]]$gamma.hat, acw_result$est[[2L]]$gamma.hat)
    
  # subject-level selective borrowing framework
  Z_mat <- diag(3000)
    
  nu.vector <- c(1.0, 2.0)
    
  # tuning the parameters (omega and lambda)
  cv_lasso_vec <- sapply(nu.vector, 
                         .cvLASSO,
                         gamma.hat = gamma_hat, 
                         bias = bias_ec, 
                         n.rct = 100, 
                         bias.score = bias_ec_score,
                         min.lambda = NULL, simplify = FALSE)

  nu_opt <- ifelse(cv_lasso_vec[[1L]]$min < cv_lasso_vec[[2L]]$min, 1L, 2L)
    
  expected <- list(cv_lasso_vec[[nu_opt]]$bias[1:1000],
                   cv_lasso_vec[[nu_opt]]$bias[1001:3000])
  
  expect_equal(.estimateLASSO(acw_result, 100, c(1.0, 2.0), NULL), expected)
})

test_that("`.estimateLASSO()` returns expected results single EC", {
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
  
  acw_result <- .estimateACW(data.rct, data.ec, aipw_result, bias)
  
  bias_ec <- list(
    acw_result$est[[1L]]$bias
  )

  n_ei <- c(1000)
  K <- 1L
  
  bias_ec_score <- list(
    acw_result$est[[1L]]$bias.score
  )
  
  gamma_hat <- c(acw_result$est[[1L]]$gamma.hat)
  
  # subject-level selective borrowing framework
  Z_mat <- diag(1000)
  
  nu.vector <- c(1.0, 2.0)
  
  # tuning the parameters (omega and lambda)
  cv_lasso_vec <- sapply(nu.vector, 
                         .cvLASSO,
                         gamma.hat = gamma_hat, 
                         bias = bias_ec, 
                         n.rct = 100, 
                         bias.score = bias_ec_score,
                         min.lambda = NULL, simplify = FALSE)
  
  nu_opt <- ifelse(cv_lasso_vec[[1L]]$min < cv_lasso_vec[[2L]]$min, 1L, 2L)
  
  expected <- list(cv_lasso_vec[[nu_opt]]$bias[1:1000])
  
  expect_equal(.estimateLASSO(acw_result, 100, c(1.0, 2.0), NULL), expected)
})