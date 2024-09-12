#' @noRd
#'
#' @param acw.result A list result. The value object returned by estimateACW.
#' @param n.rct A scalar integer. The number of participants in the RCT dataset.
#' 
#' @returns A list object. Each element the estimated bias for the participants
#'   of a single EC dataset.
#'
#' @include cvLasso.R
#' @importFrom glmnet glmnet predict.glmnet
.estimateLASSO <- function(acw.result, n.rct, nu.vector, min.lambda) {
  
  bias_ec <- lapply(acw.result$est, "[[", "bias")
  
  n_ei <- lengths(bias_ec)
  K <- length(bias_ec)
  
  bias_ec_score <- lapply(acw.result$est, "[[", "bias.score")
  
  gamma_hat <- lapply(acw.result$est, "[[", "gamma.hat") |> unlist()
  
  # subject-level selective borrowing framework
  Z_mat <- diag(sum(n_ei))
  
  if (is.null(nu.vector)) nu.vector <- c(1.0, 2.0)

  # tuning the parameters (omega and lambda)
  cv_lasso_vec <- sapply(nu.vector, 
                         .cvLASSO,
                         gamma.hat = gamma_hat, 
                         bias = bias_ec, 
                         n.rct = n.rct, 
                         bias.score = bias_ec_score,
                         min.lambda = min.lambda,
                         simplify = FALSE)
  
  trade_off <- lapply(cv_lasso_vec, "[[", "min") |> unlist()
  
  # select the one that minimize the cv error
  nu_opt <- nu.vector[which.min(trade_off)]
  
  split(cv_lasso_vec[[nu_opt]]$bias, 
        f = factor(rep(1L:K, times = n_ei))) |> unname()
  
}