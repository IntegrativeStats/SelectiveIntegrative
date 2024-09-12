#' @noRd
#'
#' @param bias.hat A numeric vector. The glmnet predicted bias.
#' @param bias A list object. The bias for all EC datasets.
#' @param bias.score A list object. The score of the bias
#' @param n.rct A scalar integer. The number of participants in the RCT dataset.
#'
#' @returns A scalar numeric.
#' 
.tradeOff <- function(bias.hat, bias, bias.score, n.rct) {
  
  K <- length(bias)
  n_ei <- lengths(bias)
  
  group_indicator <- rep(1:K, times = n_ei)
  
  zero_x <- {abs(bias.hat[-1L]) < 1e-8} |> drop()
  tau_acw_lin <- unlist(bias) * zero_x

  # for each EC, sum the non-zero bias
  bias_high <- tapply(tau_acw_lin, group_indicator,
                      function(x) { sum(x) / n.rct })
  bias_high <- sum(bias_high^2)
  
  # zeros out the elements of bias_ec_score that have zero predicted bias
  tau_acw_score_all <- unlist(lapply(bias.score, function(x) { x[-(1L:n.rct)] }))
  tau_acw_lin <- tau_acw_score_all * zero_x
  tau_acw_lin_list <- structure(split(tau_acw_lin,
                                      f = factor(group_indicator)),
                                names = NULL)
  
  tau_acw_score_recont <- mapply(function(x,y) { c(x[1L:n.rct], y) },
                                 x = bias.score, y = tau_acw_lin_list,
                                 SIMPLIFY = FALSE)

  # variance
  var_i <- unlist(mapply(function(x, y) { sum({x - mean(x)}^2) / n.rct },
                         x = tau_acw_score_recont,
                         y = n_ei, SIMPLIFY = FALSE))

  mean(var_i) + bias_high
}

#' @noRd
#' @param nu A scalar numeric. Currently hardcoded to be 1 or 2
#' @param gamma.hat A numeric vector. 
#' @param bias A list object. Each element a numeric vector.
#' @param n.rct A scalar integer. The number of participants in the RCT
#'   dataset
#' @param bias.score A list object.
#' @param min.lambda A scalar numeric.
#' 
#' @returns A numeric scalar.
#' 
#' @importFrom glmnet glmnet predict.glmnet
.cvLASSO <- function(nu, gamma.hat, bias, n.rct, bias.score, min.lambda) {
  
  # number of participants in each EC dataset
  n_ec <- lengths(bias)

  # diagonal matrix, column for each participant in the EC datasets
  Z_mat <- diag(sum(n_ec))
  
  # specify the lambda vector for cross-validation
  if (is.null(min.lambda)) min.lambda <- nrow(Z_mat)^-(1.0 + nu / 5.0) / 2.0
  max_lambda <- nrow(Z_mat)^{-0.1}
  if (max_lambda < min.lambda) max_lambda <- max(0.9, min.lambda + 0.1)

  lambda_vector <- seq(min.lambda, max_lambda, length.out = 100L)

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
  
  lasso_prediction <- predict(fit_lasso, newx = Z_mat, type = "coef")

  trade_off <-   apply(lasso_prediction, MARGIN = 2L,
                       .tradeOff,
                       bias = bias, 
                       bias.score = bias.score, 
                       n.rct = n.rct)
  
  lambda_opt_idx <- which.min(trade_off)
  
  bias_lasso <- lasso_prediction[-1L, lambda_opt_idx] |> drop()
  
  list("min" = trade_off[lambda_opt_idx],
       "bias" = bias_lasso)
}