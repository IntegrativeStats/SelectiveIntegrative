#' @noRd
#'
#' @param data.rct A list object containing elements X, A, Y, and optionally ps.
#' @param data.ec A list object. Each element contains a list of X and Y from 
#'   an external control (EC).
#' @param bias A list object. Each element containing the bias estimate for
#'   the corresponding EC dataset.
#'
#' @returns A list object containing
#' \itemize{
#'   \item est A list object. Each element the value object returned by .estimateEC
#'     for a single EC dataset.
#'   \item tau.i A list of numeric vectors. Each element of length n.rct + n.ec
#'   \item tau.hat A numeric vector of length n.ec
#'   \item tau.score A list of numeric vectors. Each element of length n.rct + n.ec
#'   \item sd.hat A numeric vector of length n.ec
#' }
#'
#' @include estimateEC.R estimateAIPW.R
#' @keywords internal
.estimateACW <- function(data.rct,
                         data.ec,
                         aipw.result,
                         bias = NULL) {
  
  if (is.null(bias)) bias <- vector(length(data.ec), mode = "list")
  
  # obtain the ACW estimator for each external control dataset
  ec_result <- mapply(.estimateEC,
                      data.ec = data.ec, bias = bias,
                      MoreArgs = list(data.rct = data.rct),
                      SIMPLIFY = FALSE)
  
  # obtain the number participants in each external control dataset
  n_ei <- lapply(data.ec, function(x) { nrow(x$X) }) |> unlist()
  n_rct <- nrow(data.rct$X)
  
  # compute the estimated standard errors for mu1_ACW-mu0_ACW
  sd_hat <- numeric(length(data.ec))
  tau_score <- list()
  tau <- list()
  tau_hat <- numeric(length(data.ec))
  for (i in seq_along(data.ec)) {
    mu <- c(aipw.result$mu$A1, rep(0.0, n_ei[i]))
    tau[[i]] <- mu - ec_result[[i]]$mu0.i
    tau_score[[i]] <- mu - ec_result[[i]]$mu0.score
    tau_hat[i] <- sum(tau[[i]], na.rm = TRUE) / n_rct
    
    mean_tau_score <- sum(tau_score[[i]], na.rm = TRUE) / n_rct
    sd_hat[i] <- {sum({tau_score[[i]] - 
                       c(rep(mean_tau_score, n_rct), rep(0.0, n_ei[i]))}^2, 
                      na.rm = TRUE) / n_rct} |> sqrt()
  }
  
  list("est" = ec_result, 
       "tau.i" = tau,
       "tau.hat" = tau_hat,
       "tau.score" = tau_score,
       "sd.hat" = sd_hat)
}