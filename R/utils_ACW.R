.var <- function(x) {
  mean({x - mean(x)}^2)
}

# calibration equation
.calEqs <- function(par, X.rct, X.ec) {

  weight_cal <- pmin(exp(par[1L] + X.ec %*% par[-1L]), 1e8) |> drop()

  left_side <- colSums(X.ec * weight_cal)
  right_side <- colSums(X.rct)
  left_side - right_side
}

.calEqsGMM <- function(par, X.rct, X.ec) {
  res <- .calEqs(par = par, X.rct = X.rct, X.ec = X.ec)
  sum(res^2)
}
