test_that("`estimateAIPW()` returns expected errors", {
  
  expect_error(.estimateAIPW(),
               "`data.rct` must be a list containing element {A, Y, ps, Y.hat.A0, Y.hat.A1}",
               fixed = TRUE)
  expect_error(.estimateAIPW(data.frame("X" = 1, "A" = 1, "Y" = 1)),
               "`data.rct` must be a list containing element {A, Y, ps, Y.hat.A0, Y.hat.A1}",
               fixed = TRUE)
  data.rct <- list()
  expect_error(.estimateAIPW(data.rct),
               "`data.rct` must be a list containing element {A, Y, ps, Y.hat.A0, Y.hat.A1}",
               fixed = TRUE)
  data.rct <- list("A" = rep(1, 10))
  expect_error(.estimateAIPW(data.rct),
               "`data.rct` must be a list containing element {A, Y, ps, Y.hat.A0, Y.hat.A1}",
               fixed = TRUE)
  data.rct$Y <- rep(1, 10)
  expect_error(.estimateAIPW(data.rct),
               "`data.rct` must be a list containing element {A, Y, ps, Y.hat.A0, Y.hat.A1}",
               fixed = TRUE)
  data.rct$ps <- rep(1, 10)
  expect_error(.estimateAIPW(data.rct),
               "`data.rct` must be a list containing element {A, Y, ps, Y.hat.A0, Y.hat.A1}",
               fixed = TRUE)
  data.rct$Y.hat.A0 <- rep(1, 10)
  expect_error(.estimateAIPW(data.rct),
               "`data.rct` must be a list containing element {A, Y, ps, Y.hat.A0, Y.hat.A1}",
               fixed = TRUE)
  data.rct$Y.hat.A1 <- rep(1, 10)

  for (i in 1L:length(data.rct)) {
    data.rct2 <- data.rct
    data.rct2[[i]] <- as.character(data.rct2[[i]])
    expect_error(.estimateAIPW(data.rct2),
                 "Elements A, Y, ps, Y.hat.A0, and Y.hat.A1 must all be numeric vectors of the same length")
  }
    
  for (i in 1L:length(data.rct)) {
    data.rct2 <- data.rct
    data.rct2[[i]] <- data.rct2[[i]][-1L]
    expect_error(.estimateAIPW(data.rct2),
                 "Elements A, Y, ps, Y.hat.A0, and Y.hat.A1 must all be numeric vectors of the same length")
  }
  
})

test_that("`estimateAIPW()` returns expected results", {
  data.rct <- withr::with_seed(
    1234, 
    list(
      "Y" = rnorm(100, 1, 8),
      "A" = rbinom(100, 1, 0.4),
      "ps" = runif(100),
      "Y.hat.A0" = rnorm(100, 1.4, 8),
      "Y.hat.A1" = rnorm(100, 1.8, 7)
    )
  )
  
  mu1 <- {data.rct$A * data.rct$Y + {data.rct$ps - data.rct$A} * data.rct$Y.hat.A1} / 
      data.rct$ps
    
  mu0 <- {{1.0 - data.rct$A} * data.rct$Y + {data.rct$A - data.rct$ps} * data.rct$Y.hat.A0} / 
      {1.0 - data.rct$ps}
    
  expected <- list("tau.hat" = sum(mu1 - mu0) / 100.0,
                   "sd.hat" = .var(mu1 - mu0) |> sqrt(),
                   "mu" = list("A0" = mu0, "A1" = mu1))
  expected$CI <- list("lower" = expected$tau.hat - 1.96 * expected$sd.hat / 10,
                      "upper" = expected$tau.hat + 1.96 * expected$sd.hat / 10)
  
  expect_equal(.estimateAIPW(data.rct), expected)

})