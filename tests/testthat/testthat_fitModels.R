test_that("`.newVariableName()` returns expected results", {
  use <- withr::with_seed(1234, paste(sample(letters, 10L, TRUE), collapse = ""))
  expect_equal(withr::with_seed(1234, .newVariableName("Y", "Y")),
               use)
  
  use <- withr::with_seed(1234, paste(sample(letters, 8L, TRUE), collapse = ""))
  expect_equal(withr::with_seed(1234, .newVariableName("Y", "Y", 8L)),
               use)

  expect_equal(withr::with_seed(1234, .newVariableName("Y2", "Y", 8L)),
               "Y")
})

test_that("`.fitEC()` returns expected errors", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2", "X3"),
         "psName" = c("X1", "X2"),
         "contName" = c("X3")))
  
  expect_error(.fitEC(data.rct, data.ec = NULL, 
                      ec.trControl = caret::trainControl(method = 'cv', number = 10L),
                      method = "gbm"),
               "`data.ec` does not appear to be in the proper format")
  
  data.ec <- withr::with_seed(
    2345,
    list("X" = matrix(rnorm(400), ncol = 2L, dimnames = list(NULL, c("X1", "X2"))),
         "Y2" = rnorm(200)))
  
  expect_error(.fitEC(data.rct, data.ec, 
                      ec.trControl = caret::trainControl(method = 'cv', number = 10L),
                      method = "gbm"),
               "`data.ec` does not appear to be in the proper format")
  
  data.ec <- withr::with_seed(
    2345,
    list("X2" = matrix(rnorm(400), ncol = 2L, dimnames = list(NULL, c("X1", "X2"))),
         "Y" = rnorm(200)))
  
  expect_error(.fitEC(data.rct, data.ec, 
                      ec.trControl = caret::trainControl(method = 'cv', number = 10L),
                      method = "gbm"),
               "`data.ec` does not appear to be in the proper format")
  
  data.ec <- withr::with_seed(
    2345,
    list("X" = matrix(rnorm(400), ncol = 2L, dimnames = list(NULL, c("X1", "X2"))),
         "Y" = rnorm(200)))
  
  expect_error(.fitEC(data.rct, data.ec, 
                      ec.trControl = caret::trainControl(method = 'cv', number = 10L),
                      method = "gbm"),
               "`data.ec` does not contain all required covariates")
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = c("X1", "X2", "X3"),
         "contName" = c("X3")))
  
  expect_error(.fitEC(data.rct, data.ec, 
                      ec.trControl = caret::trainControl(method = 'cv', number = 10L),
                      method = "gbm"),
               "`data.ec` does not contain all required covariates")
  
})

test_that("`.fitEC()` returns expected results multiple covariates", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = c("X1", "X2"),
         "contName" = c("X3"),
         "response.name" = "Y",
         "tx.name" = "A"))
  data.rct$fit.ps <- stats::glm(A ~ X1 + X2, 
                                data = data.frame("A" = data.rct$A, data.rct$X),
                                family = "binomial")
  data.rct$fit.Y <- caret::train(Y~X1+X2 + A*X3, data.frame("Y" = data.rct$Y, "A" = data.rct$A, data.rct$X),
                                 trControl = caret::trainControl(method = 'cv', number = 10L),
                                 method = "lm")
  
  
  data.ec <- withr::with_seed(
    2345,
    list("X" = matrix(rnorm(400), ncol = 2L, dimnames = list(NULL, c("X1", "X2"))),
         "Y" = rnorm(200),
         "A" = rep(0L, 200L),
         "response.name" = "Y",
         "tx.name" = "A", 
         "mainName" = c("X1", "X2"),
         "psName" = c("X1", "X2"),
         "contName" = c("X3")))
  data.ec$X <- cbind(data.ec$X, "X3" = rep(0.0, 200L))
  data.ec$ps <- stats::predict(data.rct$fit.ps, data.frame(data.ec$X), 
                               type = "response") |> drop() |> unname()

  fit_Y_ec <- withr::with_seed(
    345,
    caret::train(Y~X1+X2, data.frame("Y" = data.ec$Y, data.ec$X),
                 trControl = caret::trainControl(method = 'cv', number = 10L),
                 method = "lm"))
  
  data.ec$Y.hat <- list(
    "rct" = predict(data.rct$fit.Y, newdata = data.frame("A" = data.ec$A, data.ec$X)) |> drop() |> unname(),
    "ec" = predict(fit_Y_ec, newdata = data.frame(data.ec$X)) |> drop() |> unname()
    )

  expect_equal(withr::with_seed(345, .fitEC(data.rct, data.ec, caret::trainControl(method = 'cv', number = 10L), "lm")),
               data.ec)
  
  
  expect_no_warning(suppressMessages(.fitEC(data.rct, data.ec, caret::trainControl(method = 'cv', number = 10L), "gbm", verbose = FALSE)))
                    
  
})

test_that("`.fitEC()` returns expected results intercept only models", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = NULL,
         "psName" = NULL,
         "contName" = NULL,
         "response.name" = "Y",
         "tx.name" = "A"))
  data.rct$fit.ps <- stats::glm(A ~ 1, 
                                data = data.frame("A" = data.rct$A, data.rct$X),
                                family = "binomial")
  data.rct$fit.Y <- suppressWarnings(
    caret::train(Y~A, data.frame("Y" = data.rct$Y, "A" = data.rct$A, data.rct$X),
                 trControl = caret::trainControl(method = 'cv', number = 10L),
                 method = "lm")
  )
  
  
  data.ec <- withr::with_seed(
    2345,
    list("X" = matrix(rnorm(400), ncol = 2L, dimnames = list(NULL, c("X1", "X2"))),
         "Y" = rnorm(200),
         "response.name" = "Y",
         "tx.name" = "A", 
         "mainName" = NULL,
         "psName" = NULL,
         "contName" = NULL,
         "A" = rep(0L, 200L)))
  data.ec$ps <- stats::predict(data.rct$fit.ps, data.frame(data.ec$X), 
                               type = "response") |> drop() |> unname()
  
  fit_Y_ec <- withr::with_seed(
    345,
    suppressWarnings(
    caret::train(Y~1 + A, data.frame("Y" = data.ec$Y, "A" = data.ec$A, data.ec$X),
                 trControl = caret::trainControl(method = 'cv', number = 10L),
                 method = "lm")))

  data.ec$Y.hat <- list(
    "rct" = predict(data.rct$fit.Y, newdata = data.frame("A" = data.ec$A, data.ec$X)) |> drop() |> unname(),
    "ec" = predict(fit_Y_ec, newdata = data.frame(data.ec$X, "A" = data.ec$A)) |> drop() |> unname()
  )
  
  data.ec2 <- withr::with_seed(
    2345,
    list("X" = matrix(rnorm(400), ncol = 2L, dimnames = list(NULL, c("X1", "X2"))),
         "Y" = rnorm(200)))

  # warnings generated by caret due to intercept only model
  expect_equal(withr::with_seed(
    345, suppressWarnings(
    .fitEC(data.rct, data.ec2, caret::trainControl(method = 'cv', number = 10L), "lm"))),
               data.ec)
  
})


test_that("`.fitEC()` returns expected results intercept only models without covariates", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(0, nrow = 100L, ncol = 0L),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = NULL,
         "psName" = NULL,
         "contName" = NULL,
         "response.name" = "Y",
         "tx.name" = "A"))
  data.rct$fit.ps <- stats::glm(A ~ 1, 
                                data = data.frame("A" = data.rct$A, data.rct$X),
                                family = "binomial")
  data.rct$fit.Y <- suppressWarnings(
    caret::train(Y~A, data.frame("Y" = data.rct$Y, "A" = data.rct$A, data.rct$X),
                 trControl = caret::trainControl(method = 'cv', number = 10L),
                 method = "lm")
  )
  
  
  data.ec <- withr::with_seed(
    2345,
    list("X" = matrix(0, nrow = 400L, ncol = 0L),
         "Y" = rnorm(400),
         "response.name" = "Y",
         "tx.name" = "A", 
         "mainName" = NULL,
         "psName" = NULL,
         "contName" = NULL,
         "A" = rep(0L, 400L)))
  data.ec$ps <- stats::predict(data.rct$fit.ps, data.frame(data.ec$X), 
                               type = "response") |> drop() |> unname()
  
  fit_Y_ec <- withr::with_seed(
    345,
    suppressWarnings(
      caret::train(Y~1 + A, data.frame("Y" = data.ec$Y, "A" = data.ec$A, data.ec$X),
                   trControl = caret::trainControl(method = 'cv', number = 10L),
                   method = "lm")))
  
  data.ec$Y.hat <- list(
    "rct" = predict(data.rct$fit.Y, newdata = data.frame("A" = data.ec$A, data.ec$X)) |> drop() |> unname(),
    "ec" = predict(fit_Y_ec, newdata = data.frame(data.ec$X, "A" = data.ec$A)) |> drop() |> unname()
  )
  
  data.ec2 <- withr::with_seed(
    2345,
    list("X" = matrix(0, nrow = 400L, ncol = 0L),
         "Y" = rnorm(400)))
  
  # warnings generated by caret due to intercept only model
  expect_equal(withr::with_seed(
    345, suppressWarnings(
      .fitEC(data.rct, data.ec2, caret::trainControl(method = 'cv', number = 10L), "lm"))),
    data.ec)
  
})

test_that("`.fitRCTOutcome()` returns expected results multiple covariates", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = c("X1", "X2"),
         "contName" = c("X3"),
         "response.name" = "Y",
         "tx.name" = "A"))
  data.rct2 <- data.rct
    
  form <- Y ~ X1 + X2 + A + A:X3
  df_rct <- data.frame("Y" = data.rct$Y, "A" = data.rct$A, data.rct$X)

  fit_Y_rct <- withr::with_seed(345,
                                caret::train(form,
                                             data = df_rct,
                                             trControl = caret::trainControl(),
                                             method = "lm"))
  
  data.rct$fit.Y <- fit_Y_rct
  
  df_rct$A <- 0L
  data.rct$Y.hat.A0 <- predict(fit_Y_rct, newdata = df_rct)
  df_rct$A <- 1L
  data.rct$Y.hat.A1 <- predict(fit_Y_rct, newdata = df_rct)

  test <- withr::with_seed(345, 
                           .fitRCTOutcome(data.rct2, caret::trainControl(), "lm"))
  test$fit.Y$call <- NULL
  test$fit.Y$times <- NULL
  test$fit.Y$terms <- NULL
  
  data.rct$fit.Y$call <- NULL
  data.rct$fit.Y$times <- NULL
  data.rct$fit.Y$terms <- NULL
  
  expect_equal(test, data.rct)
  
})

test_that("`.fitRCTOutcome()` returns expected results intercept only models", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = NULL,
         "psName" = c("X1", "X2"),
         "contName" = NULL,
         "response.name" = "Y",
         "tx.name" = "A"))
  data.rct2 <- data.rct
  
  form <- Y ~ A
  df_rct <- data.frame("Y" = data.rct$Y, "A" = data.rct$A, data.rct$X)
  
  fit_Y_rct <- withr::with_seed(345, caret::train(form,
                            data = df_rct,
                            trControl = caret::trainControl(method = 'cv', number = 10L),
                            method = "lm"))
  data.rct$fit.Y <- fit_Y_rct
  
  df_rct$A <- 0L
  data.rct$Y.hat.A0 <- predict(fit_Y_rct, newdata = df_rct)
  df_rct$A <- 1L
  data.rct$Y.hat.A1 <- predict(fit_Y_rct, newdata = df_rct)
  
  test <- withr::with_seed(345, 
                           .fitRCTOutcome(data.rct2, caret::trainControl(method = 'cv', number = 10L), "lm"))
  test$fit.Y$call <- NULL
  test$fit.Y$times <- NULL
  test$fit.Y$terms <- NULL
  
  data.rct$fit.Y$call <- NULL
  data.rct$fit.Y$times <- NULL
  data.rct$fit.Y$terms <- NULL
  
  expect_equal(test, data.rct)
  
})

test_that("`.fitRCTOutcome()` returns expected results no covariates", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(0.0, nrow = 100, ncol = 0L),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = NULL,
         "psName" = NULL,
         "contName" = NULL,
         "response.name" = "Y",
         "tx.name" = "A"))
  data.rct2 <- data.rct
  
  form <- Y ~ A
  df_rct <- data.frame("Y" = data.rct$Y, "A" = data.rct$A, data.rct$X)
  
  fit_Y_rct <- withr::with_seed(345, 
                                caret::train(form,
                                             data = df_rct,
                                             trControl = caret::trainControl(),
                                             method = "lm"))
  data.rct$fit.Y <- fit_Y_rct
  
  df_rct$A <- 0L
  data.rct$Y.hat.A0 <- predict(fit_Y_rct, newdata = df_rct)
  df_rct$A <- 1L
  data.rct$Y.hat.A1 <- predict(fit_Y_rct, newdata = df_rct)
  
  test <- withr::with_seed(345, 
                           .fitRCTOutcome(data.rct2, caret::trainControl(), "lm"))
  test$fit.Y$call <- NULL
  test$fit.Y$times <- NULL
  test$fit.Y$terms <- NULL
  
  data.rct$fit.Y$call <- NULL
  data.rct$fit.Y$times <- NULL
  data.rct$fit.Y$terms <- NULL
  
  expect_equal(test, data.rct)
  
})

test_that("`.fitRCTPropensity()` returns expected results multiple covariates", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = c("X1", "X2"),
         "contName" = c("X3"),
         "response.name" = "Y",
         "tx.name" = "A"))
  
  df_rct <- data.frame(data.rct$X, "Y" = data.rct$Y, "A" = data.rct$A)
  form <- A ~ X1 + X2
  fit_ps_rct <- stats::glm(formula = form, data = df_rct, family = "binomial")
  data.rct$ps <- predict(fit_ps_rct, type = 'response') |> drop() |> unname()
  data.rct$fit.ps <- fit_ps_rct
  
  data.rct2 <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = c("X1", "X2"),
         "contName" = c("X3"),
         "response.name" = "Y",
         "tx.name" = "A"))
  
  expect_equal(.fitRCTPropensity(data.rct2, NULL),
               data.rct)
  
  data.rct$ps <- rep(0.25, 100)
  
  expect_equal(.fitRCTPropensity(data.rct2, rep(0.25, 100)), data.rct)

})

test_that("`.fitRCTPropensity()` returns expected results intercept only model", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = NULL,
         "contName" = c("X3"),
         "response.name" = "Y",
         "tx.name" = "A"))
  
  df_rct <- data.frame(data.rct$X, "Y" = data.rct$Y, "A" = data.rct$A)
  form <- A ~ 1
  fit_ps_rct <- stats::glm(formula = form, data = df_rct, family = "binomial")
  data.rct$ps <- predict(fit_ps_rct, type = 'response') |> drop() |> unname()
  data.rct$fit.ps <- fit_ps_rct
  
  data.rct2 <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = NULL,
         "contName" = c("X3"),
         "response.name" = "Y",
         "tx.name" = "A"))
  
  expect_equal(.fitRCTPropensity(data.rct2, NULL),
               data.rct)
  
  data.rct$ps <- rep(0.25, 100)
  
  expect_equal(.fitRCTPropensity(data.rct2, rep(0.25, 100)), data.rct)
  
})

test_that("`.fitRCTPropensity()` returns expected results no covariates", {
  
  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(0.0, nrow = 100, ncol = 0L),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = NULL,
         "contName" = c("X3"),
         "response.name" = "Y",
         "tx.name" = "A"))
  
  df_rct <- data.frame("Y" = data.rct$Y, "A" = data.rct$A, data.rct$X)
  form <- A ~ 1
  fit_ps_rct <- stats::glm(formula = form, data = df_rct, family = "binomial")
  data.rct$ps <- predict(fit_ps_rct, type = 'response') |> drop() |> unname()
  data.rct$fit.ps <- fit_ps_rct
  
  data.rct2 <- withr::with_seed(
    1234,
    list("X" = matrix(0.0, nrow = 100, ncol = 0L),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = NULL,
         "contName" = c("X3"),
         "response.name" = "Y",
         "tx.name" = "A"))
  
  expect_equal(.fitRCTPropensity(data.rct2, NULL),
               data.rct)
  
  data.rct$ps <- rep(0.25, 100)
  
  expect_equal(.fitRCTPropensity(data.rct2, rep(0.25, 100)), data.rct)
  
})

test_that("`.fitModels()` returns expected errors", {
  expect_error(.fitModels(), 
               "`data.rct` must be a list containing element {A, Y, X, mainName, contName, psName}",
               fixed = TRUE)
  data.rct <- withr::with_seed(
    1234, 
    list("X" = matrix(rnorm(300), nrow = 100, ncol = 3L, 
                      dimnames = list(NULL, c("X1", "X2", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("X1", "X2"),
         "psName" = 1L))
  expect_error(.fitModels(data.rct), 
               "`data.rct` must be a named list containing elements X, Y, A, mainName, contName, psName",
               fixed = TRUE)
  data.rct$contName <- "X3"
  
  expect_error(.fitModels(data.rct, data.ec = matrix(100, 10, 1)),
               "`data.ec` must be a list, each element containing a list of X and Y",
               fixed = TRUE)
  
  expect_error(.fitModels(data.rct, data.ec = list()),
               "`data.ec` must be a list, each element containing a list of X and Y",
               fixed = TRUE)
  
  expect_error(.fitModels(data.rct, data.ec = list("X" = matrix(100, 100, 1, dimnames = list(NULL, "X1")), 
                                                   "Y" = numeric(100))),
               "`each element of data.ec` must be a named list containing elements X, Y",
               fixed = TRUE)
  data.ec <- list(list("X" = matrix(100, 100, 1L, dimnames = list(NULL, "x1")),
                       "Y" = numeric(100)))
  
  expect_error(.fitModels(data.rct, data.ec),
               "`ps.rct` must be NULL or a numeric vector of length = nrow(data.rct$X)",
               fixed = TRUE)

  expect_error(.fitModels(data.rct, data.ec, ps.rct = rep(0.2, 10)),
               "`ps.rct` must be NULL or a numeric vector of length = nrow(data.rct$X)",
               fixed = TRUE)

  expect_error(.fitModels(data.rct, data.ec, ps.rct = rep(NA, 100)),
               "`ps.rct` must be NULL or a numeric vector of length = nrow(data.rct$X)",
               fixed = TRUE)
  ps.rct <- rep(0.5, 100L)
  
  expect_error(.fitModels(data.rct, data.ec, ps.rct),
               "`rct.trControl` must be a list")
  expect_error(.fitModels(data.rct, data.ec, ps.rct, rct.trControls = c("cv" = 10)),
               "`rct.trControl` must be a list")
  
  expect_error(.fitModels(data.rct, data.ec, ps.rct, list()),
               "`ec.trControl` must be a list")
  expect_error(.fitModels(data.rct, data.ec, ps.rct, list(), ec.trControls = c("cv" = 10)),
               "`ec.trControl` must be a list")
  
  expect_error(.fitModels(data.rct, data.ec, ps.rct, list(), list()),
               "`method` must be a character")
  expect_error(.fitModels(data.rct, data.ec, ps.rct, list(), list(), c("A", "B")),
               "`method` must be a character")
  expect_error(.fitModels(data.rct, data.ec, ps.rct, list(), list(), 1.0),
               "`method` must be a character")
  expect_error(.fitModels(data.rct, data.ec, ps.rct, list(), list(), list("A")),
               "`method` must be a character")
  
})

test_that("`.fitModel()` returns expected results", {

  data.rct <- withr::with_seed(
    1234,
    list("X" = matrix(rnorm(300), ncol = 3L, dimnames = list(NULL, c("Y", "A", "X3"))),
         "Y" = rnorm(100),
         "A" = rbinom(100, 1, 0.4), 
         "mainName" = c("Y", "A", "X3"),
         "psName" = c("Y", "A"),
         "contName" = c("X3")))
  
  data.ec <- withr::with_seed(
    2345,
    list("X" = matrix(rnorm(600), ncol = 3L, dimnames = list(NULL, c("Y", "A", "X3"))),
         "Y" = rnorm(200)))
  
  data.rct2 <- data.rct
  data.ec2 <- list(data.ec, data.ec)
  
  withr::with_seed(
    123, {
      data.rct$response.name <- .newVariableName(colnames(data.rct$X), "Y")
      data.rct$tx.name <- .newVariableName(c(colnames(data.rct$X), data.rct$response.name), "A")

      data.rct <- .fitRCTOutcome (data.rct = data.rct, 
                                  rct.trControl = caret::trainControl(), 
                                  method = "lm")
      
      data.rct <- .fitRCTPropensity(data.rct = data.rct, 
                                    ps.rct = NULL)

      tmp <- .fitEC(data.rct = data.rct, data.ec = data.ec,
                        ec.trControl = caret::trainControl(), method = "lm")
  
      tmp <- .fitEC(data.rct = data.rct, data.ec = data.ec,
                    ec.trControl = caret::trainControl(), method = "lm")
      
      data.ec <- list(tmp, tmp)
      data.rct$fit.Y$call <- NULL
      data.rct$fit.Y$times <- NULL
      data.rct$fit.Y$terms <- NULL
    })
  tst <- withr::with_seed(123,
                          .fitModels(data.rct2, data.ec2, NULL, caret::trainControl(), caret::trainControl(), "lm"))
  data.rct$fit.Y <- NULL
  data.rct$fit.ps <- NULL

  expect_equal(tst,
               list("data.rct" = data.rct, "data.ec" = data.ec))
  
})