#' Test if provided object is a numeric vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a numeric vector of appropriate length.
#' @keywords internal
.isNumericVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.numeric(x) && is.vector(x, mode = "numeric") && length(x) > 0L
  } else {
    is.numeric(x) && is.vector(x, mode = "numeric") && length(x) == length
  }
}

#' Test if provided object is a numeric vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#' @param nms A character vector or NULL. If NULL, no naming convention is
#'   tested; if character, all nms must be in names(x) but not necessarily
#'   in the same ordering
#'
#' @returns A logical. TRUE if x is a numeric vector of appropriate length
#'  and requested naming convention.
#' @keywords internal
.isNamedNumericVector <- function(x, length = NULL, nms = NULL) {
  tst <- is.numeric(x) && is.vector(x, mode = "numeric") && length(x) > 0L &&
    !is.null(names(x))
  if (!is.null(length)) tst <- tst && length(x) == length
  if (!is.null(nms)) tst <- tst && all(nms %in% names(x))
  tst
}

#' Test if provided object is a character vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a character vector of appropriate length.
#' @keywords internal
.isCharacterVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.character(x) && is.vector(x, mode = "character") && length(x) > 0L
  } else {
    is.character(x) && is.vector(x, mode = "character") && length(x) == length
  }
}

#' Test if provided object is an integer vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is an integer vector of appropriate length.
#' @keywords internal
.isIntegerVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.integer(x) && is.vector(x, mode = "integer") && length(x) > 0L
  } else {
    is.integer(x) && is.vector(x, mode = "integer") && length(x) == length
  }
}

#' Test if provided object is a logical vector of specified length
#'
#' @noRd
#' @param x An R object.
#' @param length An integer or NULL. If NULL, length is tested to be non-zero;
#'   if positive, length of x must match.
#'
#' @returns A logical. TRUE if x is a logical vector of appropriate length.
#' @keywords internal
.isLogicalVector <- function(x, length = NULL) {
  if (is.null(length)) {
    is.logical(x) && is.vector(x, mode = "logical") && length(x) > 0L
  } else {
    is.logical(x) && is.vector(x, mode = "logical") && length(x) == length
  }
}

#' Test if provided object is named numeric matrix
#'
#' @noRd
#' @param x An R object.
#' @param nms A character vector or NULL. If NULL, no naming convention is
#'   tested; if character, all nms must be in colnames(x) but not necessarily
#'   in the same ordering
#'
#' @returns A logical. TRUE if x is a named numeric matrix
#' @keywords internal
.isNamedNumericMatrix <- function(x, nms = NULL) {
  tst <- is.matrix(x) && !is.null(colnames(x)) && is.numeric(x) &&
    all(nchar(colnames(x)) > 0L)
  if(tst && !is.null(nms)) {
    tst <- tst & all(nms %in% colnames(x))
  }
  tst
}

#' Test if provided object is numeric matrix
#'
#' @noRd
#' @param x An R object.
#'
#' @returns A logical. TRUE if x is a numeric matrix
#' @keywords internal
.isNumericMatrix <- function(x) {
  is.matrix(x) && is.numeric(x)
}
