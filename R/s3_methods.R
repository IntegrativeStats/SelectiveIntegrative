#' @describeIn srEC Print summary of analysis
#' @param x An object of S3 class \code{SREC}.
#' @param ... Ignored
#'
#' @export
print.SREC <- function(x, ...) {
  
  print(format(summary(x), digits = 3L))
  
  invisible(x)
}

#' @describeIn srEC Summary of analysis.
#' @param object An object of S3 class \code{SREC}.
#'
#' @export
summary.SREC <- function(object, ...) {
  object$subset.idx <- NULL
  do.call(rbind, object) |> data.frame()
}