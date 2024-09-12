#' @noRd
#'
#' @param n.rct A scalar integer, the number of participants in the RCT dataset.
#' @param aipw.result A list object. The value object returned by .estimateAIPW()
#' @param acw.lasso.result A list object. The value object returned by estimateACWLASSO()
#' 
#' @returns A list object containing element tau.hat, the K+1 estimated tau
#'   and sd.hat, the K+1 estimated standard error
#'   
#' @importFrom MASS ginv
.estimateFinal <- function(n_rct, aipw.result, acw.lasso.result) {
  
  tau_score <- acw.lasso.result$tau.score # mu_AIPW - EC_score
  n_rct <- length(aipw.result$mu$A0)
  
  rct <- seq_len(n_rct)
  K <- length(tau_score)
  
  # pull out the first rct elements from each vector of tau_score
  # the cov(AIPW, ECi) terms
  tau_score_rct <- lapply(tau_score, "[", rct)
  
  term1 <- lapply(tau_score_rct, 
                  function(x) {x - sum(x, na.rm = TRUE) / length(x) })
  
  covar_acw_lasso <- lapply(term1, 
                            function(x) {
                              lapply(term1, 
                                     function(y, x) {
                                       sum(x * y, na.rm = TRUE) / n_rct
                                     }, x = x) |> unlist()
                            }) |> unlist() |> matrix(nrow = K, ncol = K)

  covar <- suppressWarnings(matrix(covar_acw_lasso, nrow = K + 1L, ncol = K + 1L))
  diag(covar) <- c(aipw.result$sd.hat^2, acw.lasso.result$sd.hat^2)

  inv_sigma <- tryCatch(MASS::ginv(covar),
                        error = function(e) {
                          stop("unable to invert Sigma matrix\n\t",
                               e$message, call. = FALSE)
                        })

  d_A_n <- rowSums(inv_sigma) / sum(inv_sigma)

  tau_final <- {c(aipw.result$tau.hat, acw.lasso.result$tau.hat) %*% d_A_n} |> drop()
  var_final <- {c(d_A_n) %*% covar %*% d_A_n} |> drop()
  sd_final <- sqrt(var_final)
  
  list("tau.hat" = tau_final,
       "sd.hat" = sd_final)
}