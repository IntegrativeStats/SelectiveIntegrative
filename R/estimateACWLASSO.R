#' @noRd
#'
#' @param data.rct A list object. The data for the RCT dataset as defined by 
#'   dataInput
#' @param data.ec A list object. Each element containing the X and Y for a single
#'   EC dataset
#' @param bias.lasso A list object. Each element a numeric vector giving the
#'   estimated bias for all participants in the EC dataset
#'   
#' @returns A list object containing tau.hat, the estimated tau for each EC
#'   dataset; sd.hat, the estimated sd for each EC dataset; and tau.score
#'   A list of numeric vectors. Each element of length n.rct + n.ec.
#'
#' @include estimateACW.R
.estimateACWLASSO <- function(data.rct, data.ec, bias.lasso) {
  
  acw_result <- .estimateACW(data.rct = data.rct,
                             data.ec = data.ec,
                             bias = bias.lasso)
  n_rct <- nrow(data.rct$X)
  
  tau_hat <- mapply(function(x, y) {
                      zero_bias <- as.numeric(abs(y) < 1e-8)
                      sum(x$tau * c(rep(1, n_rct), zero_bias), na.rm = TRUE) / n_rct
                    },
                    x = acw_result$tau.i,
                    y = bias.lasso, SIMPLIFY = FALSE) |> unlist()

  list("tau.hat" = tau_hat,
       "sd.hat" = sqrt(acw_result$sd.hat^2),
       "tau.score" = acw_result$tau.score)
}