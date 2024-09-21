.newVariableName <- function(existing, proposed, max.length = 10L) {
  variable_name <- proposed
  while (variable_name %in% existing) {
    variable_name <- paste(sample(letters, max.length, TRUE), collapse = "")
  }
  variable_name
}

.fitEC <- function(data.rct, data.ec, ec.trControl, method, ...) {
  
  if (!is.vector(data.ec, "list") ||
      !all(c("X", "Y") %in% names(data.ec))) {
    stop("`data.ec` does not appear to be in the proper format", call. = FALSE)
  }
  
  if (!all(c(data.rct$mainName, data.rct$psName) %in% c(colnames(data.ec$X)))) {
    stop("`data.ec` does not contain all required covariates", call. = FALSE)
  }
  
  # Add zero columns if any contName variables are missing
  if (!all(data.rct$contName %in% colnames(data.ec$X))) {
    are_missing <- !{data.rct$contName %in% colnames(data.ec$X)}
    data.ec$X <- cbind(matrix(0.0, nrow(data.ec$X), sum(are_missing),
                              dimnames = list(NULL, data.rct$contName[are_missing])),
                       data.ec$X)
  }
  
  df_ec <- data.frame(data.ec$X)
  df_ec[[data.rct$tx.name]] <- data.ec$A

  # will just be easier if everything has the same info...
  carry_over <- c("response.name", "tx.name", "mainName", "psName", "contName")
  data.ec[carry_over] <- data.rct[carry_over]
  
  data.ec$A <- rep(0L, nrow(data.ec$X))
  
  data.ec$ps <- predict(data.rct$fit.ps, newdata = df_ec, 
                        type = "response") |> drop() |> unname()
  
  # Note that A will be all zeros for ec data, including here to
  # ensure that caret can handle intercept only models. The NA
  # values of the coefficients do not cause issue
  form <- paste(data.rct$response.name, "~", 
                ifelse(is.null(data.rct$mainName), "XXXX1111",
                       paste(data.rct$mainName, collapse = " + "))) |> stats::as.formula()
  
  df_ec <- data.frame(data.ec$X)
  df_ec[[data.rct$response.name]] <- data.ec$Y
  df_ec[[data.rct$tx.name]] <- data.ec$A

  if (is.null(data.rct$mainName)) {
    # Have to fool caret sometimes for the intercept only models
    df_ec$XXXX1111 <- 1.0
    fit_Y_ec <- tryCatch(caret::train(form,
                                      data = df_ec,
                                      trControl = ec.trControl,
                                      method = method, intercept = FALSE, ...),
                          error = function(e) {
                            stop("error encountered in caret::train() of EC outcome\n\t",
                                 e$message, call. = FALSE)
                          })
    
  } else {
    fit_Y_ec <- tryCatch(caret::train(form,
                                      data = df_ec,
                                      trControl = ec.trControl,
                                      method = method, ...),
                         error = function(e) {
                           stop("error encountered in caret::train() of EC outcome\n\t",
                                e$message, call. = FALSE)
                         })
  }
  
  data.ec$Y.hat <- list(
    "rct" = predict(data.rct$fit.Y, newdata = df_ec) |> drop() |> unname(),
    "ec" = predict(fit_Y_ec, newdata = df_ec) |> drop() |> unname()
  )
  
  data.ec
}

.fitRCTOutcome <- function(data.rct, rct.trControl, method, ...) {
  
  # create model and fit the RCT data
  form <- paste(data.rct$response.name, "~", 
                ifelse(is.null(data.rct$mainName), "",
                       paste(paste(data.rct$mainName, collapse = " + "),
                             "+")), 
                data.rct$tx.name, 
                ifelse(is.null(data.rct$contName), "",
                       paste("+",
                             data.rct$tx.name, ":(", 
                             paste(data.rct$contName, collapse = " + "),
                             ")"))) |> stats::as.formula()
  
  df_rct <- data.frame(data.rct$X)
  df_rct[[data.rct$response.name]] <- data.rct$Y
  df_rct[[data.rct$tx.name]] <- data.rct$A
  
  data.rct$fit.Y <- tryCatch(caret::train(form,
                                          data = df_rct,
                                          trControl = rct.trControl,
                                          method = method, ...),
                        error = function(e) {
                          stop("error encountered in caret::train() of RCT outcome\n\t",
                               e$message, call. = FALSE)
                        })
  
  df_rct[[data.rct$tx.name]] <- 0L
  data.rct$Y.hat.A0 <- predict(data.rct$fit.Y, newdata = df_rct)
  df_rct[[data.rct$tx.name]] <- 1L
  data.rct$Y.hat.A1 <- predict(data.rct$fit.Y, newdata = df_rct)

  data.rct
}

.fitRCTPropensity <- function(data.rct, ps.rct) {
 
  df_rct <- data.frame(data.rct$X)
  df_rct[[data.rct$response.name]] <- data.rct$Y
  df_rct[[data.rct$tx.name]] <- data.rct$A

  form <- paste(data.rct$tx.name, "~", 
                ifelse(is.null(data.rct$psName), "1", 
                       paste(data.rct$psName, collapse = " + "))) |> stats::as.formula()
  fit_ps_rct <- tryCatch(stats::glm(formula = form, data = df_rct, family = "binomial"),
                         error = function(e) {
                           stop("error encountered in glm fit of RCT propensity\n\t",
                                e$message, call. = FALSE)
                         })
  
  # if propensity not specified by user, set based on model predictions
  if (is.null(ps.rct)) {
    data.rct$ps <- predict(fit_ps_rct, type = 'response') |> drop() |> unname()
  } else {
    data.rct$ps <- ps.rct
  }
  
  data.rct$fit.ps <- fit_ps_rct
  data.rct
  
}

#' @noRd
#' 
#' @param data.rct A list object containing elements X, A, Y, and optionally ps.
#' @param data.ec NULL or a list object. If a list, each element contains the X 
#'   and Y from a single external control (EC).
#' @param ps.rct NULL or a numeric vector. Optional input providing a vector of
#'   known propensity scores P(A=1|X) for the RCT dataset. If not provided,
#'   it will be estimated using the model defined in \code{data.rct$psName}.
#' @param rct.trControl A list object. The control list for caret::train for
#'   the RCT dataset
#' @param ec.trControl A list object. The control list for caret::train for
#'   the EC datasets
#'
#' @returns data.rct and data.ec augmented with fitted models
#'
#' @import caret
#' @importFrom stats as.formula glm predict
#' @include dataInput.R stopTests.R
.fitModels <- function(data.rct, data.ec, ps.rct, 
                       rct.trControl, ec.trControl, method, ...) {

  stopifnot(
    "`data.rct` must be a list containing element {A, Y, X, mainName, contName, psName}" =
      !missing(data.rct) && .isDI(data.rct, "data.rct"),
    "`data.ec` must be a list, each element containing a list of X and Y" =
      !missing(data.ec) && {is.null(data.ec) || 
      {is.vector(data.ec, "list") && length(data.ec) >= 1L &&
          all(lapply(data.ec, .isReducedDI, "each element of data.ec") |> unlist())}},
    "`ps.rct` must be NULL or a numeric vector of length = nrow(data.rct$X)" =
      !missing(ps.rct) &&
      {is.null(ps.rct) || .isNumericVector(ps.rct, nrow(data.rct$X))},
    "`rct.trControl` must be a list" = !missing(rct.trControl) &&
      is.vector(rct.trControl, "list"),
    "`ec.trControl` must be a list" = !missing(ec.trControl) &&
      is.vector(ec.trControl, "list"),
    "`method` must be a character" = !missing(method) && .isCharacterVector(method, 1L)
  )

  # identify a unique response variable name
  response_name <- .newVariableName(colnames(data.rct$X), "Y")
  # identify a unique treatment variable name
  tx_name <- .newVariableName(c(colnames(data.rct$X), response_name), "A")

  data.rct$response.name <- response_name
  data.rct$tx.name <- tx_name
  
  data.rct <- .fitRCTOutcome (data.rct = data.rct, 
                              rct.trControl = rct.trControl, 
                              method = method, ...)
  
  # estimate propensity score (if EC provided or ps not provided)
  if(is.null(ps.rct) || !is.null(data.ec)) {
    
    data.rct <- .fitRCTPropensity(data.rct = data.rct, 
                                  ps.rct = ps.rct)
    
    # allow for the possibility of providing a list of the X and Y when
    # only 1 EC provided
    if (length(data.ec) >= 2L && all(c("X", "Y") %in% names(data.ec))) {
      data.ec <- list(data.ec)
    }
    
    for (i in seq_along(data.ec)) {
      data.ec[[i]] <- .fitEC(data.rct = data.rct, data.ec = data.ec[[i]],
                             ec.trControl = ec.trControl, method = method, ...)
    }
    
  }
  data.rct$fit.Y <- NULL
  data.rct$fit.ps <- NULL
  
  
  list("data.rct" = data.rct, "data.ec" = data.ec)
}