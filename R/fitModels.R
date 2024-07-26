.newVariableName <- function(existing, proposed, max.length = 10L) {
  variable_name <- proposed
  while (variable_name %in% existing) {
    variable_name <- sample(letters, max.length, TRUE)
  }
  variable_name
}

#' @noRd
#' 
#' @param data.rct A list object containing elements X, A, Y, and optionally ps.
#' @param data.ec A list object. Each element contains a list of X and Y from 
#'   an external control (EC).
#' @param rct.trControl A list object. The control list for caret::train for
#'   the RCT dataset
#' @param ec.trControl A list object. The control list for caret::train for
#'   the EC datasets
#'
#' @returns data.rct and data.ec augmented with fitted models
#'
#' @import caret
#' @importFrom stats as.formula glm predict
.fitModels <- function(data.rct, data.ec, rct.trControl, ec.trControl) {
 
  stopifnot(
    "`data.rct` must be a list containing element {A, Y, X, mainName, contName, psName}" =
      !missing(data.rct) && .isDI(data.rct),
    "`data.ec` must be a list, each element containing a list of X and Y" =
      !missing(data.ec) && all(lapply(data.ec, .isReducedDI) |> unlist()),
    "`rct.trControl` must be a list" = !missing(rct.trControl) &&
      is.vector(rct.trControl, "list"),
    "`ec.trControl` must be a list" = !missing(ec.trControl) &&
      is.vector(ec.trControl, "list")
  )
  
  # identify a unique response variable name
  response_name <- .newVariableName(colnames(data.rct$X), "Y")
  # identify a unique treatment variable name
  tx_name <- .newVariableName(c(colnames(data.rct$X), response_name), "A")

  # create model and fit the RCT data
  form <- paste(response_name, "~", 
                paste(data.rct$mainName, collapse = " + "),
                "+", tx_name, "+",
                tx_name, ":(", 
                paste(data.rct$contName, collapse = " + "),
                ")") |> stats::as.formula()
  
  df_rct <- data.frame("Y" = data.rct$Y, "A" = data.rct$A, data.rct$X)
  colnames(df_rct)[1L:2L] <- c(response_name, tx_name)
  
  data.rct$response.name <- response_name
  data.rct$tx.name <- tx_name
  
  fit_Y_rct <- tryCatch(caret::train(form,
                                     data = df_rct,
                                     trControl = rct.trControl,
                                     method = method, ...),
                        error = function(e) {
                          stop("error encountered in caret::train() of RCT outcome\n\t",
                               e$message, call. = FALSE)
                        })
  
  df_rct[[tx_name]] <- 0L
  data.rct$Y.hat.A0 <- predict(fit_Y_rct, newdata = df_rct)
  df_rct[[tx_name]] <- 1L
  data.rct$Y.hat.A1 <- predict(fit_Y_rct, newdata = df_rct)
  df_rct[[tx_name]] <- data.rct$A
  
  
  # estimate propensity score (if EC provided or ps not provided)
  if(is.null(data.rct$ps) || !is.null(data.ec)) {
    
    form <- paste(tx_name, "~", 
                  paste(data.rct$psName, collapse = " + ")) |> stats::as.formula()
    fit_ps_rct <- tryCatch(stats::glm(formula = form, data = df_rct, family = binomial()),
                           error = function(e) {
                             stop("error encountered in glm fit of RCT propensity\n\t",
                                  e$message, call. = FALSE)
                           })
    
    # if propensity not specified by user, set based on model predictions
    if (is.null(data.rct$ps)) data.rct$ps <- predict(fit_ps_rct, type = 'response') |> drop() |> unname()
    
    for (i in seq_along(data.ec)) {
      
      if (!all(c(data.rct$mainName, data.rct$psName) %in% colnames(data.ec[[i]]$X))) {
        stop("`data.ec` does not contain all required covariates", call. = FALSE)
      }
      
      # Add zero columns if any contName variables are missing
      if (!all(data.rct$contName %in% colnames(data.ec[[i]]$X))) {
        which_are_missing <- !{data.rct$contName %in% colnames(data.ec[[i]]$X)}
        data.ec[[i]]$X <- cbind(matrix(0.0, nrow(data.ec[[i]]$X), sum(which_are_missing),
                                       dimnames = list(NULL, data.rct$contName[which_are_missing])),
                                data.ec[[i]]$X)
      }
      
      df_ec <- data.frame(A = 0L, data.ec[[i]]$X)
      colnames(df_ec)[1L] <- tx_name
      
      # will just be easier if everything has the same info...
      data.ec[[i]]$response.name <- response_name
      data.ec[[i]]$tx.name <- tx_name
      data.ec[[i]]$mainName <- data.rct$mainName
      data.ec[[i]]$psName <- data.rct$psName
      data.ec[[i]]$contName <- data.rct$contName
      
      data.ec[[i]]$ps <- predict(fit_ps_rct, newdata = df_ec, type = "response") |> drop() |> unname()
      
      form <- paste(data.ec$response.name, "~", 
                    paste(data.ec$mainName, collapse = " + ")) |> stats::as.formula()
      
      df_ec <- data.frame("Y" = data.ec$Y, "A" = 0L, data.ec$X)
      colnames(df_ec)[1L:2L] <- c(data.ec$response.name, data.ec$tx.name)
      
      fit_Y_ec <- tryCatch(caret::train(form,
                                        data = df_ec,
                                        trControl = ec.trControl,
                                        method = method, ...),
                           error = function(e) {
                             stop("error encountered in caret::train() of EC outcome\n\t",
                                  e$message, call. = FALSE)
                           })
      
      data.ec[[i]]$Y.hat <- list("rct" = predict(fit_Y_rct, newdata = df_ec) |> drop() |> unname(),
                                 "ec" = predict(fit_Y_ec, newdata = df_ec) |> drop() |> unname())
    }
    
  }
  
  list(data.rct = data.rct, data.rwe = data.rwe)
}