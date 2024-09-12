#' @noRd
#' @param data.rct A list object containing {X, mainName}
#' @param data.ec A list object containing {X, mainName}
#' 
#' @returns A number vector of length(mainName) + 1
#' 
#' @include utils_ACW.R
#' @importFrom stats optim
#' @importFrom nleqslv nleqslv
#' 
#' @keywords internal
.lambdaHat <- function(data.rct, data.ec) {
  
  stopifnot(
    "`data.rct` must be a list containing {X, mainName}" = 
      !missing(data.rct) && is.vector(data.rct, "list") &&
      all(c("X", "mainName") %in% names(data.rct)),
    "`data.ec` must be a list containing {X, mainName}" = 
      !missing(data.ec) && is.vector(data.ec, "list") &&
      all(c("X", "mainName") %in% names(data.ec))
  )

  p_main_effect <- length(data.rct$mainName) + 1L
  
  n_ec <- nrow(data.ec$X)
  
  if (p_main_effect < n_ec) {
    lambda_hat <- tryCatch(
      nleqslv::nleqslv(x = rep(0.0, p_main_effect), # initial points
                       fn = .calEqs,
                       jac = .calEqsGradient,
                       X.rct = data.rct$X[, data.rct$mainName, drop = FALSE], 
                       X.ec = data.ec$X[, data.ec$mainName, drop = FALSE], 
                       control = list(trace = FALSE,
                                      maxit = 500L),
                       method = 'Newton'),
      warning = function(w) {
        stats::optim(par = rep(0.0, p_main_effect),
                     fn = .calEqsGMM,
                     X.rct = data.rct$X[, data.rct$mainName, drop = FALSE], 
                     X.ec = data.ec$X[, data.ec$mainName, drop = FALSE], 
                     control = list(trace = TRUE,
                                    abstol = 1e-8,
                                    maxit = 500L),
                     method = 'BFGS')
      },
      error = function(e) {
        stop("error encountered in nleqslv::nleqslv\n\t",
             e$message, call. = FALSE)
      })
    
    if ("par" %in% names(lambda_hat)) {
      if (lambda_hat$convergence != 0) {
        warning("stats::optim returned a non-0 convergence flag: ", 
                lambda_hat$convergence, call. = FALSE)
      }
      
      lambda_hat <- lambda_hat$par
    } else {
      if (lambda_hat$termcd != 1) {
        warning("nleqslv::nleqslv returned termination code: ",
                lambda_hat$termcd, call. = FALSE)
      }
      lambda_hat <- lambda_hat$x
    }
    
  } else {
    lambda_hat <-  tryCatch(stats::optim(par = rep(0.0, p_main_effect),
                                         fn = .calEqsGMM,
                                         X.rct = data.rct$X[, data.rct$mainName, drop = FALSE], 
                                         X.ec = data.ec$X[, data.ec$mainName, drop = FALSE], 
                                         control = list(trace = TRUE,
                                                        abstol = 1e-8,
                                                        maxit = 500L),
                                         method = 'BFGS'),
                            error = function(e) {
                              stop("error encountered in stats::optim\n\t",
                                   e$message, call. = FALSE)
                            })
    
    if (lambda_hat$convergence != 0) {
      warning("stats::optim returned a non-0 convergence flag: ", 
              lambda_hat$convergence, call. = FALSE)
    }
    
    lambda_hat <- lambda_hat$par
  }
  lambda_hat
}

.SACW <- function(X.rct, X.ec, q.hat.ec) {
  # {n_rct + n_ec x p_me}
  S_acw_q <- rbind(-X.rct, X.ec * q.hat.ec)
  # {p_me x p_me}
  dot_S_acw_q <- t(q.hat.ec * X.ec) %*% X.ec / nrow(X.rct)
  
  # {p_me x p_me}
  inv_inf <- tryCatch(MASS::ginv(dot_S_acw_q),
                      error = function(e) {
                        stop("error encountered in evaluating inverse of influence function\n\t",
                             e$message, call. = FALSE)
                      })
  
  list("S" = S_acw_q, "inv.inf" = inv_inf)
  
}

#' @noRd
#' @param data.rct A list object containing elements {A, Y, X, mainName, ps, Y.hat.A0}
#' @param data.ec A list object containing {X, Y} for a single EC dataset
#' @param lambda.hat A numeric vector of length p_main_effect
#' @param r.X A scalar numeric.
#' @param zero.bias A logical vector of length data.ec$Y
#' 
#' @returns A list containing numeric vectors q.hat.ec (length n.ec), 
#'   mu0.i (length n.rct + n.ec), and mu0.score (length n.rct + n.ec).
#' 
#' @importFrom MASS ginv
#'
# NOTE: IT IS ASSUMED THAT lambda.hat IS IN THE ORDER OF 
# Intercept, mainName
.influenceFunction <- function(data.rct, data.ec, lambda.hat, r.X, zero.bias) {
  
  stopifnot(
    "`data.rct` must be a list containing {X, A, Y, mainName, ps, Y.hat.A0}" = 
      !missing(data.rct) && is.vector(data.rct, "list") && 
      all(c("X", "A", "Y", "mainName", "ps", "Y.hat.A0") %in% names(data.rct)),
    "`data.ec` must be a list containing {X, Y, mainName, ps, Y.hat}" = !missing(data.ec) &&
      is.vector(data.ec, "list") && all(c("X", "Y", "mainName", "ps", "Y.hat") %in% names(data.ec)),
    "`data.ec$Y.hat` must be a list containing element {rct}" =
      is.vector(data.ec$Y.hat, "list") && all(c("rct") %in% names(data.ec$Y.hat)),
    "`lambda.hat` must be a numeric vector of lenght p_main_effects" = 
      !missing(lambda.hat) && .isNumericVector(lambda.hat, length(data.rct$mainName) + 1L),
    "`r.X` must be a positive scalar numeric" = !missing(r.X) &&
      .isNumericVector(r.X, 1L) && r.X >= 0.0,
    "`zero.bias` must be a logical vector" = !missing(zero.bias) &&
      .isLogicalVector(zero.bias, length(data.ec$Y))
  )
  
  n_rct <- nrow(data.rct$X)
  
  # {n_rct x p_me}
  ME_rct_with_intercept <- cbind(1.0, data.rct$X[, data.rct$mainName, drop = FALSE])
  # {n_ec x p_me}
  ME_ec_with_intercept <- cbind(1.0, data.ec$X[, data.ec$mainName, drop = FALSE])
  
  # {n_ec}
  q_hat_ec <- pmin(exp(ME_ec_with_intercept %*% lambda.hat), 1e8) |> drop()
  # {n_rct}
  q_hat_rct <- pmin(exp(ME_rct_with_intercept %*% lambda.hat), 1e8) |> drop()

  S_results <- .SACW(X.rct = ME_rct_with_intercept, 
                     X.ec = ME_ec_with_intercept, 
                     q.hat.ec = q_hat_ec)
  
  deno_rct <- r.X + {1.0 - data.rct$ps} * q_hat_rct
  numo_rct <- q_hat_rct *  {1.0 - data.rct$A} * {data.rct$Y - data.rct$Y.hat.A0}
  tmp <- r.X * numo_rct / deno_rct^2 * ME_rct_with_intercept
  dot_mu0_rct_q <- colSums(tmp, na.rm = TRUE) / n_rct
  
  deno_ec <- r.X * zero.bias + {1.0 - data.ec$ps} * q_hat_ec
  numo_ec <- q_hat_ec * r.X * zero.bias * {data.ec$Y - data.ec$Y.hat$rct}
  tmp <- r.X * numo_ec / deno_ec^2 * ME_ec_with_intercept
  dot_mu0_ec_q <- colSums(tmp, na.rm = TRUE) / n_rct
  
  # {n_rct + n_ec}
  mu0_i <- c(data.rct$Y.hat.A0 + numo_rct / deno_rct, numo_ec / deno_ec)
  
  # {n_rct + n_ec}
  mu0_score <-  mu0_i - 
    drop(tcrossprod(S_results$S, {dot_mu0_rct_q + dot_mu0_ec_q} %*% S_results$inv.inf))
  
  list("q.hat.ec" = q_hat_ec,
       "mu0.i" = mu0_i,
       "mu0.score" = mu0_score)
  
}

#' @noRd
#' @param data.rct A list object containing elements X, A, Y, and optionally ps.
#' @param data.ec A list object containing elemenst X and Y for a single EC
#' @param bias An optional numeric vector. If given, length must match the
#'   dimensions dictated by data.ec.
#'   
#' @returns A list containing numeric vectors 
#'   q.hat.ec (length n.ec), 
#'   mu0.i (length n.rct + n.ec), 
#'   mu0.score (length n.rct + n.ec),
#'   bias (length n.ec),
#'   bias.score (length n.rct + n.ec),
#    gamma.hat (length n.ec), and
#'   r.X (length 1)
#'   
#' @include utils_ACW.R
#' @keywords internal
.estimateEC <- function(data.rct, data.ec, bias = NULL) {
  
  stopifnot(
    "`data.rct` must be a list containing elements {X, A, Y, mainName, ps, Y.hat.A0}" =
      !missing(data.rct) && is.vector(data.rct, "list") &&
      all (c("A", "X", "Y", "mainName", "ps", "Y.hat.A0") %in% names(data.rct)),
    "`data.ec` must be a list containing elements {X, Y, mainName, ps, Y.hat}" =
      !missing(data.ec) && is.vector(data.ec, "list") &&
      all (c("X", "Y", "mainName", "ps", "Y.hat") %in% names(data.ec)),
    "`data.ec$Y.hat` must be a list containing elements {rct, ec}" =
      is.vector(data.ec$Y.hat, "list") && all(c("ec", "rct") %in% names(data.ec$Y.hat)),
    "`bias` must be NULL or a numeric vector" = is.null(bias) || 
      .isNumericVector(bias, length(data.ec$Y))
  )
  
  n_rct <- nrow(data.rct$X)
  n_ec <- nrow(data.ec$X)
  
  # vector {p_main_effects}
  lambda_hat <- .lambdaHat(data.rct = data.rct, data.ec = data.ec)

  # the default value for bias
  if (is.null(bias)) { 
    bias <- rep(0.0, n_ec)
    zero_bias <- rep(TRUE, n_ec)
  } else {
    zero_bias <- abs(bias) < 1e-8
  }
  
  tmp <- if (any(zero_bias)) {
    mean({data.ec$Y.hat$ec[zero_bias] - mean(data.ec$Y[zero_bias])}^2)
  } else {
    1e6 
  }
  
  r_X <- mean({data.rct$Y.hat.A0[data.rct$A == 0L] - mean(data.rct$Y[data.rct$A == 0L])}^2) / tmp

  # construct the pseudo-observations for the bias parameters
  # {n_ec}
  bias <- data.ec$Y - data.ec$Y.hat$rct

  # compute the score for estimating q
  # {n_ec + n_rct}
  bias_i_score <- c(rep(0.0, n_rct), bias)
  
  # construct the influence function
  mu0 <- .influenceFunction(data.rct = data.rct, 
                            data.ec = data.ec, 
                            lambda.hat = lambda_hat, 
                            r.X = r_X, 
                            zero.bias = zero_bias)
  
  c(mu0,
    list("bias" = bias,
         "bias.score" = bias_i_score,
         "gamma.hat" = data.ec$Y.hat$ec - data.ec$Y.hat$rct, 
         "r.X" = r_X))
}