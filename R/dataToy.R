#' Toy Continuous Outcome Dataset
#'
#' These datasets are provided only to facilitate examples. They are not based
#'   on or representative of any real-world applications.
#'
#' @name selectiveToy
#' @rdname selectiveToy
#' @aliases selectiveToy.rct selectiveToy.rwe
#'
#' @usage data("selectiveToy")
#'
#' @format selectiveToy provides two datasets. The selectiveToy.rct
#'   2986 participant records; selectiveToy.rwe 1014 participant records. Each
#'   data.frame provides the following:
#' \itemize{
#'   \item Y: A continuous outcome.
#'   \item X1: A continuous covariate.
#'   \item X2: A continuous covariate.
#'   \item A: A binary treatment variable
#' }
#'
#' @keywords datasets
NULL

#' Generate continuous data for example
#'
#' Code used to create package datasets. Provided for the convenience of future
#'   developer, not intended for use by users. Default settings are those
#'   used to generate the data provided with the package. The code is not
#'   robustly tested.
#'
#' @noRd
#' @param seed n An integer. The random seed for data generation
#' @param n.rct An integer. The size of the clinical trial dataset.
#' @param n.rwe An integer. The size of the real-world evidence dataset.
#' @param omega A numeric. The unmeasured confounding for EC.
#'
#' @returns A list containing trivial datasets provided with the package
#'
#' @importFrom stats rbinom rnorm uniroot
#' @keywords internal
.generateToyContData <- function(seed = 2333L, n.rct = 3000L, n.rwe = 1000L, omega = 3.0) {
  
  set.seed(seed)
  
  N <- n.rct + n.rwe
  X1 <- stats::rnorm(N)
  X2 <- stats::rnorm(N)
  
  esp <- function(alpha0, X, beta) {
    exp(alpha0 + X %*% beta) / (1.0 + exp(alpha0 + X %*% beta)) |> drop()
  }
  
  f <- function(alpha0, X, N, n, beta) {
    mean(esp(alpha0, X, beta)) * N - n
  }
  
  alpha0.opt <- stats::uniroot(f, 
                               interval = c(-50, 10),
                               beta = c(-2.0, -2.0),
                               X = cbind(X1, X2), N = N, n = n.rct)$root

  eS <- esp(alpha0.opt, cbind(X1, X2), beta = c(-2.0, -2.0))
  print(alpha0.opt)
  print(summary(eS))
  delta <- stats::rbinom(N, size = 1, prob = eS)
  
  X.rt <- cbind(X1, X2)[delta == 1L, ]

  ## generate the treatment assignment with marginal probability P.A
  P.A <- 2.0 / 3.0
  eta0.opt <- stats::uniroot(f, interval = c(-50, 10),
                             beta = c(-1.0, -1.0, -1.0),
                             X = cbind(1.0, X.rt),
                             N = 1.0,
                             n = P.A)$root
  eA <- esp(eta0.opt,cbind(1.0, X.rt), beta = c(-1.0, -1.0, -1.0))
  A.rt <- stats::rbinom(nrow(X.rt), size = 1, prob = eA)
  
  ## generate the observed outcomes for RT
  Y.rt <- {1.0 +  1.0 + X.rt %*% c(1, 1) + A.rt * X.rt %*% c(.3, .3) + 
    stats::rnorm(nrow(X.rt)) + 
    omega * stats::rnorm(nrow(X.rt))} |> drop()

  selectiveToy.rct <- data.frame("X1" = X.rt[, 1L], 
                                 "X2" = X.rt[, 2L],
                                 "Y" = Y.rt, 
                                 "A" = A.rt)
  
  # generate the external control population
  X.ec <- cbind(X1, X2)[delta == 0L, ]
  ## generate the observed outcomes for EC (possibly confounded)
  Y.ec <- {1.0 +  1.0 + X.ec %*% c(1, 1) + 
      stats::rnorm(nrow(X.ec)) +
      omega * stats::rnorm(nrow(X.ec), mean = 1)} |> drop()

  selectiveToy.rwe <- data.frame("X1" = X.ec[, 1L], 
                                 "X2" = X.ec[, 2L],
                                 "Y" = Y.ec,
                                 "A" = 0L)

  list("rct" = selectiveToy.rct, "rwe" = selectiveToy.rwe)
}
