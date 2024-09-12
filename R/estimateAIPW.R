#' Treatment Effect Estimated Using AIPW
#' 
#' @noRd
#' 
#' @param data.rct A list object containing elements A, Y, and model
#'   predictions Y.hat.A0, Y.hat.A1, and ps for the RCT dataset.
#'
#' @returns A list object with elements
#' \itemize{
#'   \item tau.hat: The estimated treatment effect.
#'   \item sd.hat: The standard error of the estimated treatment effect
#'   \item mu: A list containing mu for A = 1 and A = 0
#'   \item CI: A list giving the "lower" and "upper" confidence intervals
#' }
#' 
#' @include utils_ACW.R
.estimateAIPW <- function(data.rct) {
  
  stopifnot(
    "`data.rct` must be a list containing element {A, Y, ps, Y.hat.A0, Y.hat.A1}" =
      !missing(data.rct) && 
      is.vector(data.rct, mode = "list") &&
      all(c("A", "Y", "ps", "Y.hat.A0", "Y.hat.A1") %in% names(data.rct)),
    "Elements A, Y, ps, Y.hat.A0, and Y.hat.A1 must all be numeric vectors of the same length" =
      all(lapply(data.rct[c("A", "Y", "ps", "Y.hat.A0", "Y.hat.A1")], is.vector, mode = "numeric") |> unlist()) &&
      all(lengths(data.rct[c("A", "Y", "ps", "Y.hat.A0", "Y.hat.A1")]) == length(data.rct$A))
  )

  n_rct <- length(data.rct$Y)
  
  mu1 <- {data.rct$A * data.rct$Y + {data.rct$ps - data.rct$A} * data.rct$Y.hat.A1} / 
    data.rct$ps
  
  mu0 <- {{1.0 - data.rct$A} * data.rct$Y + {data.rct$A - data.rct$ps} * data.rct$Y.hat.A0} / 
    {1.0 - data.rct$ps}
  
  tau_i <-  mu1 - mu0
  tau_hat <- mean(tau_i)
  sd_hat <- .var(tau_i) |> sqrt()
  
  tau_upper <- tau_hat + 1.96 * sd_hat / sqrt(n_rct)
  tau_lower <- tau_hat - 1.96 * sd_hat / sqrt(n_rct)
  
  list("tau.hat" = tau_hat,
       "sd.hat" = sd_hat,
       "mu" = list("A0" = mu0, "A1" = mu1),
       "CI" = list("lower" = tau_lower,
                   "upper" = tau_upper))
    
}