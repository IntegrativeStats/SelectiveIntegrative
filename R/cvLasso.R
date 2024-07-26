#' @noRd
#'
#' @param bias.hat A numeric vector. The glmnet predicted bias.
#' @param bias A list object. The bias for all EC datasets.
#' @param bias.score A list object. The score of the bias
#' @param n.rct A scalar integer. The number of participants in the RCT dataset.
#'
#' @returns A scalar numeric.
#' 
.tradeOff2 <- function(bias.hat, bias, bias.score, n.rct) {
  
  # the number of EC datasets
  K <- length(bias)
  
  # the number of participants in each EC dataset
  n_ec <- lengths(bias)
  
  # vector to resplit bias from vector -> list
  group_indicator <- rep(1L:K, times = n_ec)
  
  # remove bias values that are predicted to be zero
  zero_bias <- abs(bias.hat[-1L]) < 1e-8
  tau_acw_lin <- unlist(bias) * zero_bias
  
  # STH I expected this to be n_ec or the number of nonzero n_ec
  # the average bias in each EC dataset
  bias_high <- tapply(tau_acw_lin, group_indicator, function(x) { sum(x) / n.rct })
  # STH I expected this to be K not n.rct
  bias_high <- sum(bias_high^2) / n.rct
  
  # remove the RCT participants from each bias.score
  tau_acw_score_all <- unlist(lapply(bias.score, function(x) { x[-(1L:n.rct)] }))
  tau_acw_lin <- tau_acw_score_all * zero_bias
  # STH I expected this to be n_ec or the number of nonzero n_ec
  var_high <- tapply(tau_acw_lin, group_indicator,
                     function(x) { sum({x - sum(x) / n.rct}^2) / n.rct })
            
  # keep only the RCT participants from each bias.score
  tau_acw_score_rct <- lapply(bias.score, function(x) { x[(1L:n.rct)] })
  var_ci <- lapply(tau_acw_score_rct, .var) |> unlist()
            
  bias_high + mean(var_high) + mean(var_ci)
}

#' @noRd
#' @param nu A scalar numeric. Currently hardcoded to be 1 or 2
#' @param gamma.hat A numeric vector. 
#' @param bias A list object. Each element a numeric vector.
#' @param n.rct A scalar integer. The number of participants in the RCT
#'   dataset
#' @param bias.h.score A list object.
#' 
#' @returns A numeric scalar.
#' 
#' @import glmnet glmnet predict.glmnet
.cvLASSO <- function(nu, gamma.hat, bias, n.rct, bias.score) {
  
  # number of participants in each EC dataset
  n_ec <- lengths(bias)

  # diagonal matrix, column for each participant in the EC datasets
  Z_mat <- diag(sum(n_ec))
  
  # specify the lambda vector for cross-validation
  lambda_vector <- seq(nrow(Z_mat)^-(1.0 + nu / 5.0) / 2.0,
                       nrow(Z_mat)^(-0.1),
                       length.out = 100L)
  
  # construct the adaptive weights (subject-level)
  b_k_w <- c(1.0 / abs(gamma.hat)^nu)
  
  # the penalty factors are re-scaled with sum equal to `nvars`
  b_k_w <- b_k_w / sum(b_k_w) * ncol(Z_mat)
  
  # replace any NA number
  b_k_w[is.na(b_k_w)] <- Inf
  
  # fit the adaptive lasso penalized estimates for the subject-level bias
  fit_lasso <- tryCatch(glmnet::glmnet(Z_mat, unlist(bias), 
                                       family = "gaussian",
                                       alpha = 1,
                                       penalty.factor = b_k_w,
                                       intercept = FALSE,
                                       standardize = FALSE,
                                       standardize.response = FALSE,
                                       lambda = lambda_vector),
                        error = function(e) {
                          stop("error encounter in glmnet for nu=", nu, "\n\t",
                               e$message, call. = FALSE)
                        })
  
  # EC dataset membership for all participants
  
  lasso_prediction <- glmnet::predict.glmnet(fit_lasso, newx = Z_mat, type = "coef")
  
  trade_off <-   apply(lasso_prediction, MARGIN = 2L,
                       .tradeOff2,
                       bias = bias, 
                       bias.score = bias.score, 
                       n.rct = n.rct)

  min(trade_off)
}