.tradeOff <- function(bias.hat, bias, bias.score, n.rct) {
  
  K <- length(bias)
  n_ei <- lengths(bias)
  
  group_indicator <- rep(1:K, times = n_ei)
  
  zero_x <- {abs(bias.hat[-1L]) < 1e-8} |> drop()
  
  tau_acw_lin <- unlist(bias) * zero_x
  
  # for each EC, sum the non-zero bias
  bias_high <- tapply(tau_acw_lin, group.indicator,
                      function(x) { sum(x) / n.rct })
  bias_high <- sum(bias_high^2)
  
  # zeros out the elements of bias_ec_score that have zero predicted bias
  tau_acw_score_all <- unlist(lapply(bias_ec_score, function(x) { x[-(1::n.rct)] }))
  tau_acw_lin <- tau_acw_score_all * zero_x
  tau_acw_lin_list <- structure(split(tau_acw_lin,
                                      f = factor(group_indicator)),
                                names = NULL)
  
  tau_acw_score_recont <- mapply(function(x,y) { c(x[1L:n.rct], y) },
                                 x = bias_ec_score, y = tau_acw_lin_list,
                                 SIMPLIFY = FALSE)
  # variance
  var_i <- unlist(mapply(function(x, y) { sum({x - mean(x)}^2) / n.rct },
                         x = tau_acw_score_recont,
                         y = est_acw_result$n_ei, SIMPLIFY = FALSE))
  mean(var_i) + bias_high
}


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
.estimateLASSO <- function(acw.result, n.rct) {
  
  bias_ec <- lapply(acw.result$est, "[[", "bias")
  
  n_ei <- lengths(bias_ec)
  K <- length(bias_ec)
  
  bias_ec_score <- lapply(acw.result$est, "[[", "bias.score")
  
  gamma_hat <- lapply(acw.result$est, "[[", "gamma.hat") |> unlist()
  
  # subject-level selective borrowing framework
  Z_mat <- diag(sum(n_ei))
  
  
  # tuning the parameters (omega and lambda)
  nu.vector <- c(1.0, 2.0)
  cv_lasso_vec <- sapply(nu.vector, 
                         .cvLASSO,
                         gamma.hat = gamma_hat, 
                         bias.h = bias_ec, 
                         K = K, 
                         n_ei = n_ei, 
                         n.rct = n.rct, 
                         bias.h.score = bias_ec_score)
  
  # select the one that minimize the cv error
  nu_opt <- nu.vector[which.min(cv_lasso_vec)]
  
  b_k_OLS <- c(1.0 / abs(gamma_hat)^nu_opt)
  b_k_w <- b_k_OLS / sum(b_k_OLS) * ncol(Z_mat)
  b_k_w[is.na(b_k_w)] <- 1e10
  
  # nu_opt/3 due to a smaller order of convergence rate of the machine learning models
  lambda_vector <- c(seq(nrow(Z_mat)^-(-.001 + nu_opt / 3.0) / 2.0, # change to 3 for conservative convergence rate
                         nrow(Z_mat)^(-0.001),
                         length.out = 100L))
  
  fit_lasso <- glmnet::glmnet(Z.mat, bias_ec, family = 'gaussian',
                              alpha = 1,
                              penalty.factor = b_k_w,
                              intercept = FALSE,
                              standardize = FALSE,
                              standardize.response = FALSE,
                              lambda = lambda_vector)
  
  lasso_prediction <- glmnet::predict(fit_lasso, newx = Z.mat, type = 'coef')
  
  trade_off <- apply(lasso_prediction, 2,
                     .tradeOff,
                     bias = bias_ec, 
                     bias.score = bias_ec_score, 
                     n.rct = n.rct)

  lambda_opt_idx <- which.min(trade_off)
  
  # STH If it turns out that trade off expression should not be different,
  # can eliminate the code above and just call the cvLasso with the new
  # lambda -- will need to add a new argument because of the nu_opt / 3
  # vs the nu / 5 in cvLasso
  
  # obtain the penalized estimates for the subject-level bias
  bias_lasso <- lasso_prediction[-1L, lambda_opt_idx] |> drop()
  
  split(bias_lasso, f = factor(rep(1L:K, times = n_ei))) |> unname()
  
}