test_that("`.estimateFinal()` returns expected results multiple EC", {
  
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
  
  acw_lasso_result <- .estimateACWLASSO(data.rct, data.ec, aipw_result, bias)
  
  tau_score <- acw_lasso_result$tau.score
  tau_score_rct <- list()
  tau_score_rct[[1L]] <- tau_score[[1L]][1:100]
  tau_score_rct[[2L]] <- tau_score[[2L]][1:100]

  term1 <- list()
  term1[[1L]] <- tau_score_rct[[1L]] - sum(tau_score_rct[[1L]], na.rm = TRUE) / 100
  term1[[2L]] <- tau_score_rct[[2L]] - sum(tau_score_rct[[2L]], na.rm = TRUE) / 100
  
  covar_acw_lasso <- matrix(0.0, 2, 2)
  covar_acw_lasso[1L, 1L] <- sum(term1[[1L]] * term1[[1L]]) / 100
  covar_acw_lasso[1L, 2L] <- sum(term1[[1L]] * term1[[2L]]) / 100
  covar_acw_lasso[2L, 1L] <- sum(term1[[2L]] * term1[[1L]]) / 100
  covar_acw_lasso[2L, 2L] <- sum(term1[[2L]] * term1[[2L]]) / 100

  sigma2_group <- matrix(0, 3, 3)
  sigma2_group[1L, 1L] <- aipw_result$sd.hat^2
  sigma2_group[1L, 2L] <- covar_acw_lasso[2L, 2L]
  sigma2_group[1L, 3L] <- covar_acw_lasso[1L, 2L]
  
  sigma2_group[2L, 1L] <- covar_acw_lasso[2L, 1L]
  sigma2_group[2L, 2L] <- acw_lasso_result$sd.hat[1L]^2
  sigma2_group[2L, 3L] <- covar_acw_lasso[2L, 2L]
  
  sigma2_group[3L, 1L] <- covar_acw_lasso[1L, 2L]
  sigma2_group[3L, 2L] <- covar_acw_lasso[2L, 1L]
  sigma2_group[3L, 3L] <- acw_lasso_result$sd.hat[2L]^2

  inv_sigma <- MASS::ginv(sigma2_group)
  d_A_n <- rowSums(inv_sigma) / sum(inv_sigma)

  tau_final <- {c(aipw_result$tau.hat, acw_lasso_result$tau.hat) %*% d_A_n} |> drop()
  var_final <- {c(d_A_n) %*% sigma2_group %*% d_A_n} |> drop()
  sd_final <- sqrt(var_final)
    
  expected <- list("tau.hat" = tau_final,
                   "sd.hat" = sd_final)

  expect_equal(.estimateFinal(100, aipw_result, acw_lasso_result), expected)
})

test_that("`.estimateFinal()` returns expected results single EC", {
  
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
  
  acw_lasso_result <- .estimateACWLASSO(data.rct, data.ec, aipw_result, bias)
  
  tau_score <- acw_lasso_result$tau.score
  tau_score_rct <- list()
  tau_score_rct[[1L]] <- tau_score[[1L]][1:100]

  term1 <- list()
  term1[[1L]] <- tau_score_rct[[1L]] - sum(tau_score_rct[[1L]], na.rm = TRUE) / 100

  covar_acw_lasso <- matrix(0.0, 1, 1)
  covar_acw_lasso[1L, 1L] <- sum(term1[[1L]] * term1[[1L]]) / 100

  sigma2_group <- matrix(0, 2, 2)
  sigma2_group[1L:1L, 2L:2L] <- covar_acw_lasso
  sigma2_group[2L:2L, 1L:1L] <- covar_acw_lasso
  sigma2_group[1L, 1L] <- aipw_result$sd.hat^2
  sigma2_group[2L, 2L] <- acw_lasso_result$sd.hat[1L]^2

  inv_sigma <- MASS::ginv(sigma2_group)
  d_A_n <- rowSums(inv_sigma) / sum(inv_sigma)
  
  tau_final <- {c(aipw_result$tau.hat, acw_lasso_result$tau.hat) %*% d_A_n} |> drop()
  var_final <- {c(d_A_n) %*% sigma2_group %*% d_A_n} |> drop()
  sd_final <- sqrt(var_final)
  
  expected <- list("tau.hat" = tau_final,
                   "sd.hat" = sd_final)
  
  expect_equal(.estimateFinal(100, aipw_result, acw_lasso_result), expected)
})