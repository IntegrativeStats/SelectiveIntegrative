#' Extract response variable from provided formula
#'
#' @noRd
#' @param formula A formula object of the form LHS ~ RHS
#' @param data A data.frame object containing at a minimum the covariates and
#'   response variable of formula
#'
#' @returns A vector object.
#'
#' @importFrom stats model.frame model.response
#' @keywords internal
.extractResponseVariable <- function(formula, data) {

  stopifnot(
    "`formula` must be a formula of the form LHS ~ RHS" = !missing(formula) &&
      inherits(formula, "formula") && length(formula) == 3L,
    "`data` must be a data.frame" = !missing(data) && is.data.frame(data))

  response <- tryCatch(stats::model.response(stats::model.frame(formula, data)),
                       error = function(e) {
                         stop("unable to extract response variable\n\t",
                              e$message, call. = FALSE)
                       })
  if (is.null(response)) {
    stop("response variable not in provided formula", format(formula),
         call. = FALSE)
  }
  response
}

#' Process `outcome.model` to identify main effects and contrast components
#'
#' @noRd
#' @param outcome.model A formula object of the form LHS ~ RHS
#' @param ps.model A formula object of the form LHS ~ RHS. Used to extract
#'   variable name of treatment.
#'
#' @returns A list object. Element ME contains a formula representation of the
#'   main effect model. Element CO contains a formula representation of the
#'   Z for the HTE.
#'
#' @importFrom stats as.formula terms
#' @keywords internal
.separateOutcomeModel <- function(outcome.model, ps.model) {

  stopifnot(
    "`outcome.model` must be a formula of the form LHS ~ RHS" = !missing(outcome.model) &&
      inherits(outcome.model, "formula") && length(outcome.model) == 3L,
    "`ps.model` must be a formula of the form LHS ~ RHS" = !missing(ps.model) &&
      inherits(ps.model, "formula") && length(ps.model) == 3L)

  avar <- ps.model[[2L]] |> as.character()

  outcome_factors <- attr(stats::terms(outcome.model), "factors")

  if (length(outcome_factors) == 0L) {
    warning(rep("! ", 10L), "\n",
            "`outcome.model` is intercept only; \n",
            "assuming intercept only for main effects model and Z for HTE\n",
            rep("! ", 10L), call. = FALSE)
    outcome_factors <- matrix(1L, nrow = 1L, ncol = 1L,
                              dimnames = list(avar, avar))

  }

  # main effects are those terms that do not interact with avar
  factors_a_row <- which(rownames(outcome_factors) == avar)
  if (length(factors_a_row) == 0L) {
    warning(rep("! ", 10L), "\n",
            "unable to identify treatment variable in `outcome.model`\n",
            "assuming intercept only Z for HTE\n",
            rep("! ", 10L), call. = FALSE)

    # add a row and column for the treatment variable
    outcome_factors <- rbind(outcome_factors,
                             integer(ncol(outcome_factors)))
    rownames(outcome_factors)[nrow(outcome_factors)] <- avar

    outcome_factors <- cbind(outcome_factors,
                             c(rep(0L, nrow(outcome_factors) - 1L), 1L))
    colnames(outcome_factors)[ncol(outcome_factors)] <- avar

    factors_a_row <- which(rownames(outcome_factors) == avar)
  } else if (length(factors_a_row) > 1L) {
    stop("multiple rows attributed to ", avar,
         " in attr(stats::terms(outcome.model), 'factor')",
         " verify input", call. = FALSE)
  }

  # separate the factors matrix into main effects and Z
  avar_related_cols <- which(outcome_factors[factors_a_row, ] > 0)
  avar_related_cols_names <- colnames(outcome_factors)[avar_related_cols]

  me_factors <- outcome_factors[, -avar_related_cols, drop = FALSE]
  co_factors <- outcome_factors[, avar_related_cols, drop = FALSE]

  # remove avar row
  me_factors <- me_factors[-factors_a_row, , drop = FALSE]
  co_factors <- co_factors[-factors_a_row, , drop = FALSE]

  co_factors <- co_factors[, colnames(co_factors) != avar, drop = FALSE]

  if (any(colSums(me_factors) == 0L) || any(colSums(co_factors) == 0L)) {
    stop("failed to separate main effects and Z; contact developer", call. = FALSE)
  }

  # column names of me_factors defines main effects model
  if (is.null(colnames(me_factors))) {
    me_model <- stats::as.formula("~ 1")
  } else  {
    me_model <- stats::as.formula(paste("~ 1",
                                        paste(colnames(me_factors), collapse = "+"),
                                        sep = "+"))
  }

  # have to be more clever for contrasts to remove avar from expression
  co_model <- "1"
  i <- 1L
  while (i <= ncol(co_factors)) {
    if (all(co_factors[, i] == 0L)) next
    co_model <- c(co_model,
                  paste(rownames(co_factors)[co_factors[, i] > 0L], collapse = ":"))
    i <- i + 1L
  }
  co_model <- stats::as.formula(paste("~", paste(co_model, collapse = "+")))

  list("ME" = me_model, "CO" = co_model)

}

#' Prepare Required Data Inputs to \code{srEC()}
#'
#' Provided a data.frame and models for the outcome and propensity regressions,
#'   generates an list of the format required by \code{srEC()}.
#'
#' @param data A data.frame object. Must contain the outcome of interest,
#'   the treatment and all model covariates.
#' @param outcome.model A formula object. Must be of the form LHS ~ RHS. The
#'   LHS is the outcome of interest. The RHS defines the main effects model and
#'   the Z covariates of the treatment effect function; i.e., of the form
#'   LHS ~ main_effects + treatment_variable * Z.
#' @param ps.model A formula object. Must be of the form LHS ~ RHS. The
#'   LHS is the treatment variable. The RHS defines the propensity score model.
#'
#' @returns A list object appropriately formatted for use in \code{elasticHTE()}.
#'   \itemize{
#'     \item \code{X}: The covariate matrix excluding an intercept. It
#'       contains all covariates required of the outcome and propensity
#'       score models.
#'     \item \code{Y} The response vector.
#'     \item \code{A} The treatment vector.
#'     \item \code{mainName} The column headers of X defining the main effects
#'       component of the outcome model or 1, indicating an intercept-only model.
#'     \item \code{contName} The column headers of X defining Z or 1, indicating
#'       an intercept-only model.
#'     \item \code{psName} The column headers of X defining the propensity score
#'       model or 1, indicating an intercept-only model.
#'   }
#'
#' @examples
#' data("selectiveToy")
#'
#' data_object <- dataInput(selectiveToy.rct,
#'                          outcome.model = Y ~ X1 + X2 * A,
#'                          ps.model = A ~ X1 * X2)
#'
#' @importFrom stats as.formula model.matrix
#' @export
dataInput <- function(data, outcome.model, ps.model) {

  if (!missing(outcome.model) && !inherits(outcome.model, "formula")) {
    outcome.model <- tryCatch(stats::as.formula(outcome.model),
                              error = function(e) {
                                stop("unable to create formula from `outcome.model`",
                                     e$message, call. = FALSE)
                              })
  }

  if (!missing(ps.model) && !inherits(ps.model, "formula")) {
    ps.model <- tryCatch(stats::as.formula(ps.model),
                         error = function(e) {
                           stop("unable to create formula from `ps.model`",
                                e$message, call. = FALSE)
                           })
  }

  stopifnot(
    "`data` must be a data.frame" = !missing(data) && is.data.frame(data),
    "`outcome.model` must be a formula of the form LHS ~ RHS" = !missing(outcome.model) &&
      inherits(outcome.model, "formula") && length(outcome.model) == 3L,
    "`ps.model` must be a formula of the form LHS ~ RHS" = !missing(ps.model) &&
      inherits(ps.model, "formula") && length(ps.model) == 3L)

  input_obj <- list()

  # extract response variable from outcome.model; store as Y
  input_obj$Y <- .extractResponseVariable(outcome.model, data)

  # extract response variable from ps.model; store as A
  # because we are requiring a 3 element formula, can safely extract the 2nd element
  input_obj$A <- .extractResponseVariable(ps.model, data)

  # process outcome.model to identify main effects and contrast components
  outcome_components <- .separateOutcomeModel(outcome.model, ps.model)

  # Create design matrix and use column headers to define mainName
  X_me <- stats::model.matrix(outcome_components$ME, data)[, -1L, drop = FALSE]
  input_obj$mainName <- colnames(X_me)
  if (length(input_obj$mainName) == 0L) input_obj$mainName <- 1L

  # Create design matrix and use column headers to define contName
  X_co <- stats::model.matrix(outcome_components$CO, data)[, -1L, drop = FALSE]
  input_obj$contName <- colnames(X_co)
  if (length(input_obj$contName) == 0L) input_obj$contName <- 1L

  # Create design matrix and use column headers to define psName

  X_ps <- stats::model.matrix(ps.model, data)
  if (attr(stats::terms(ps.model), "intercept") == 1L) {
    X_ps <- X_ps[, -1L, drop = FALSE]
  }
  input_obj$psName <- colnames(X_ps)
  if (length(input_obj$psName) == 0L) input_obj$psName <- 1L

  input_obj$X <- X_me
  input_obj$X <- cbind(input_obj$X, X_co[, !{colnames(X_co) %in% colnames(input_obj$X)}, drop = FALSE])
  input_obj$X <- cbind(input_obj$X, X_ps[, !{colnames(X_ps) %in% colnames(input_obj$X)}, drop = FALSE])

  input_obj
}

#' @include stopTests.R
.isDI <- function(object, object.name = "object") {

  if (missing(object)) stop("`object` must be provided", call. = FALSE)

  required_elements <- c("X", "Y", "A", "mainName", "contName", "psName")

  if (!is.vector(object, mode = "list") || !all(required_elements %in% names(object))) {
    stop("`", object.name, "` must be a named list containing elements ",
         paste(required_elements, collapse = ", "), call. = FALSE)
  }

  if (!{.isNamedNumericMatrix(object$X) || ncol(object$X) == 0L}) {
    stop("`", object.name, "$X` must be a matrix with column names",
         call. = FALSE)
  }

  if (!{{.isNumericVector(object$mainName, 1L) && isTRUE(all.equal(object$mainName, 1))} ||
        {.isCharacterVector(object$mainName) && all(object$mainName %in% colnames(object$X))}}) {
    stop("`", object.name, "$mainName` must be 1L or a character vector of X column headers",
         call. = FALSE)
  }

  if (!{{.isNumericVector(object$contName, 1L) && isTRUE(all.equal(object$contName, 1))} ||
        {.isCharacterVector(object$contName) && all(object$contName %in% colnames(object$X))}}) {
    stop("`", object.name, "$contName` must be 1L or a character vector of X column headers",
         call. = FALSE)
  }

  if (!{{.isNumericVector(object$psName, 1L) && isTRUE(all.equal(object$psName, 1))} ||
        {.isCharacterVector(object$psName) && all(object$psName %in% colnames(object$X))}}) {
    stop("`", object.name, "$psName` must be 1L or a character vector of X column headers",
         call. = FALSE)
  }
  TRUE
}

#' @include stopTests.R
.isDIProcessed <- function(object, object.name = "object") {
  
  if (missing(object)) stop("`object` must be provided", call. = FALSE)
  
  required_elements <- c("X", "Y", "A")
  
  if (!is.vector(object, mode = "list") || !all(required_elements %in% names(object))) {
    stop("`", object.name, "` must be a named list containing elements ",
         paste(required_elements, collapse = ", "), call. = FALSE)
  }
  
  if (!{.isNamedNumericMatrix(object$X) || ncol(object$X) == 0L}) {
    stop("`", object.name, "$X` must be a matrix with column names",
         call. = FALSE)
  }
  
  if (!{is.null(object$mainName) ||
      {.isCharacterVector(object$mainName) && all(object$mainName %in% colnames(object$X))}}) {
    stop("`", object.name, "$mainName` must be NULL or a character vector of X column headers",
         call. = FALSE)
  }
  
  if (!{is.null(object$contName) ||
      {.isCharacterVector(object$contName) && all(object$contName %in% colnames(object$X))}}) {
    stop("`", object.name, "$contName` must be NULL or a character vector of X column headers",
         call. = FALSE)
  }
  
  if (!{is.null(object$psName) ||
      {.isCharacterVector(object$psName) && all(object$psName %in% colnames(object$X))}}) {
    stop("`", object.name, "$psName` must be NULL or a character vector of X column headers",
         call. = FALSE)
  }
  TRUE
}



.isReducedDI <- function(object, object.name = "object") {

  if (missing(object)) stop("`object` must be provided", call. = FALSE)
  
  required_elements <- c("X", "Y")
  
  if (!is.vector(object, mode = "list") || !all(required_elements %in% names(object))) {
    stop("`", object.name, "` must be a named list containing elements ",
         paste(required_elements, collapse = ", "), call. = FALSE)
  }

  if (is.null(object$X) || !{.isNamedNumericMatrix(object$X) || ncol(object$X) == 0L}) {
    stop("`", object.name, "$X` must be a matrix with column names",
         call. = FALSE)
  }
  
  TRUE
}