.var <- function(x, n = length(x)) {
  sum({x - sum(x) / n}^2) / n
}

# calibration equation
#' @noRd
#' @param par A numeric vector. The current parameter estimates
#' @param X.rct A numeric matrix. The number of columns is length(par) - 1L
#' @param X.ec A numeric matrix. The number of columns is length(par) - 1L
#' @returns A numeric vector of length = length(par)
.calEqs <- function(par, X.rct, X.ec) {
  
  weight_cal <- pmin(exp(par[1L] + X.ec %*% par[-1L]), 1e8) |> drop()

  left_side <- c(sum(weight_cal), colSums(X.ec * weight_cal))
  right_side <- c(nrow(X.rct), colSums(X.rct))

  left_side - right_side
}

.calEqsGradient <- function(par, X.rct, X.ec) {
  weight_cal <- pmin(exp(par[1L] + X.ec %*% par[-1L]), 1e8) |> drop()
  t(weight_cal * cbind(1.0, X.ec)) %*% cbind(1.0, X.ec)
}

#' @noRd
#' @param par A numeric vector. The current parameter estimates
#' @param X.rct A numeric matrix. The number of columns is length(par) - 1L
#' @param X.ec A numeric matrix. The number of columns is length(par) - 1L
#' @returns A scalar numeric.
.calEqsGMM <- function(par, X.rct, X.ec) {
  res <- .calEqs(par = par, X.rct = X.rct, X.ec = X.ec)
  sum(res^2)
}

.adjustModelCoding <- function(model.names, available.covariates) {
  if (is.null(model.names)) {
    model.names <- available.covariates
  } else if (is.numeric(model.names)) {
    model.names <- NULL
  } else {
    if (!all(model.names %in% available.covariates)) {
      stop("unrecognized model covariate provided", call. = FALSE)
    }
  }
  model.names
}
