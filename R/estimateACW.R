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
.estimateACW <- function(data.rct,
                         data.ec,
                         aipw.result,
                         bias = list(NULL)) {
  
  # obtain the ACW estimator for each external control dataset
  data.ec_i <- mapply(.estimateEC,
                      data.ec = data.ec, bias.b = bias,
                      MoreArgs = list(data.rct = data.rct),
                      SIMPLIFY = FALSE)
  
  # obtain the number participants in each external control dataset
  n_ei <- lapply(data.ec, function(x) { nrow(x$X) }) |> unlist()
  
  # compute the estimated standard errors for mu1_ACW-mu0_ACW
  sd_hat <- numeric(length(data.ec))
  tau_score <- list()
  tau <- list()
  tau_hat <- numeric(length(data.ec))
  for (i in seq_along(data.ec)) {
    mu <- c(aipw.result$mu$A1, rep(0.0, n_ei[i]))
    tau[[i]] <- mu - data.ec_i[[i]]$mu0.i
    tau_score[[i]] <- mu - data.ec_i[[i]]$mu0.score
    # STH verify that this is n_rct (might be n_rct + n_ec)
    tau_hat[i] <- sum(tau[[i]], na.rm = TRUE) / n_rct
    # STH verify that this is n_rct (might be n_rct + n_ec)
    sd_hat[i] <- {sum({tau_score[[i]] - 
                       c(rep(tau.hat[i], n_rct), rep(0.0, n_ei[i]))}^2, 
                      na.rm = TRUE) / n_rct} |> sqrt()
  }
  
  list("est" = data.ec_i, 
       "tau.i" = tau,
       "tau.hat" = tau_hat,
       "tau.score" = tau_score,
       "sd.hat" = sd_hat)
}