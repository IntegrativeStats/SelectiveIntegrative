test_that("`srEC()` returns expected errors", {
  expect_error(srEC(),
               "`data.rct` must be provided")
  
  expect_error(srEC(list("X" = matrix(0, 10, 2))),
               "`data.rct` must be a named list containing elements X, Y, A, mainName, contName, psName")
  
  data.rct <- dataInput(data.frame("Y" = 1:100,
                                   "A" = 1:100,
                                   "X1" = 1:100,
                                   "X2" = 1:100), 
                        Y ~ X1*X2 + A*X2 + A*X1, A ~ X1 + X2)
  
  expect_error(srEC(data.rct = data.rct, data.ec = matrix(0, 10, 2)),
               "`data.ec` must be a list")
  
  expect_error(srEC(data.rct, ps.rct = numeric(50)),
               "`ps.rct` must be NULL or a numeric vector of length = nrow(data.rct$X)",
               fixed = TRUE)
  
  expect_error(srEC(data.rct, rct.trControl = c("maxit" = 50)),
               "`rct.trControl` must be a list()")

  expect_error(srEC(data.rct, ec.trControl = c("maxit" = 50)),
               "`ec.trControl` must be a list()")
  
  expect_error(srEC(data.rct, method = 2),
               "`method` must be a character")
  
  expect_error(srEC(data.rct, nu.vector = "2"),
               "`nu.vector` must be a numeric vector")

  expect_error(srEC(data.rct, min.lambda = NA),
               "`min.lambda` must be NULL, scalar numeric or a numeric vector")
  expect_error(srEC(data.rct, min.lambda = 1:3),
               "`min.lambda` must be NULL, scalar numeric or a numeric vector")
  expect_error(srEC(data.rct, min.lambda = c("a", "b")),
               "`min.lambda` must be NULL, scalar numeric or a numeric vector")
  
})

test_that("`srEC()` returns expected result", {

  data.rct <- dataInput(withr::with_seed(2345,
                                         data.frame("Y" = runif(100),
                                                    "A" = rbinom(100, 1, 0.3),
                                                    "X1" = rnorm(100, 1, 2),
                                                    "X2" = rnorm(100, 2, 1))), 
                        Y ~ X1*X2 + A*X2 + A*X1, A ~ X1 + X2)
  
  rct_trControl <- caret::trainControl()

  updated_data_objects <- .fitModels(data.rct = data.rct, data.ec = NULL,
                                     ps.rct = NULL,
                                     rct.trControl = rct_trControl,
                                     ec.trControl = rct_trControl,
                                     method = "lm", metric = "MAE")
  
  aipw_result <- .estimateAIPW(data.rct = updated_data_objects$data.rct)
  
  expected <- list("AIPW" = aipw_result[c("tau.hat", "sd.hat", "CI")]) 
  expect_equal(srEC(data.rct, rct.trControl = rct_trControl, 
                    method = "lm", metric = "MAE"),
               expected)
  
})


test_that("`srEC()` returns expected result", {
  
  data.rct <- dataInput(withr::with_seed(34531,
                                         data.frame("Y" = runif(100),
                                                    "A" = rbinom(100, 1, 0.3),
                                                    "X1" = rnorm(100, 1, 2),
                                                    "X2" = rnorm(100, 2, 1))), 
                        Y ~ X1*X2 + A*X2 + A*X1, A ~ X1 + X2)
  
  data.ec <- list(dataInput(withr::with_seed(3234124,
                                             data.frame("Y" = runif(1000),
                                                        "A" = rep(0, 1000),
                                                        "X1" = rnorm(1000, 1, 2),
                                                        "X2" = rnorm(1000, 2, 1))), 
                            Y ~ X1*X2 + A*X2 + A*X1, A ~ X1 + X2))
  
  rct_trControl <- caret::trainControl()
  ec_trControl <- caret::trainControl()
  
  updated_data_objects <- .fitModels(data.rct = data.rct, data.ec = data.ec,
                                     ps.rct = NULL,
                                     rct.trControl = rct_trControl,
                                     ec.trControl = ec_trControl,
                                     method = "lm", metric = "MAE")
  
  aipw_result <- .estimateAIPW(data.rct = updated_data_objects$data.rct)
  
  acw_result <- .estimateACW(data.rct = updated_data_objects$data.rct, 
                             data.ec = updated_data_objects$data.ec,
                             aipw.result = aipw_result)
  
  bias_lasso <- .estimateLASSO(acw.result = acw_result, n.rct = 100,
                               nu.vector = c(1, 2), min.lambda = NULL)
  
  # obtain the indices for borrowing
  ec_idx_lasso <- {abs(unlist(bias_lasso)) < 1e-8} |> which() |> unname()
  
  if (length(ec_idx_lasso) == 0L) {
    # if no external controls are selected
    expected <- list("aipw" = aipw_result[c("tau.hat", "sd.hat", "CI")],
                 "acw" = list("tau.hat" = acw_result$tau.hat,
                              "sd.hat" = acw_result$sd.hat),
                 "subset.idx" = integer(0L)) 
    print("HERE")
  } else {
  
  acw_lasso_result <- .estimateACWLASSO(data.rct = updated_data_objects$data.rct,
                                        data.ec = updated_data_objects$data.ec,
                                        aipw.result = aipw_result,
                                        bias.lasso = bias_lasso)
  
  acw_final <- .estimateFinal(aipw.result = aipw_result, 
                              acw.lasso.result = acw_lasso_result)
  
  expected <- list("aipw" = aipw_result[c("tau.hat", "sd.hat")],
       "acw" = acw_result[c("tau.hat", "sd.hat")],
       "acw.lasso" = acw_lasso_result[c("tau.hat", "sd.hat")],
       "acw.final" = acw_final,
       "subset.idx" = ec_idx_lasso)
  
  }
  
  expect_equal(srEC(data.rct, data.ec[[1L]], rct.trControl = rct_trControl, 
                    ec.trControl = ec_trControl, method = "lm", metric = "MAE"),
               expected)
  
})