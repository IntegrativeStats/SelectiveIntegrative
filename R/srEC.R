#' Selective and Robust External-Control Borrowing Strategy for Evaluating Treatment Effect
#'
#' @description
#'  A penalized dynamic integrative framework to augment a randomized clinical
#'  trial (RT) with one or more external control (EC) datasets, in which the 
#'  subject-level compatibility of each EC dataset is assessed by a well-crafted 
#'  penalty (e.g., adaptive lasso penalty). The parameter of interest is the 
#'  average treatment effect.
#'
#' @param data.rct The value object returned by \code{dataInput()} for the
#'   data from a randomized clinical trial (RCT). See \link{dataInput} for
#'   further details.
#' @param data.ec NULL or a list object. If NULL, only AIPW estimates will be
#'   calculated. If a list, each element is a list containing the X and Y for
#'   for a single external control (EC) dataset for which treatment A=0. 
#' @param ... Other options used to fit the predictive models. Passed on to 
#'   [caret::train()].
#' @param ps.rct NULL or a numeric vector. Optional input providing a vector of
#'   known propensity scores P(A=1|X) for the RCT dataset. If not provided,
#'   it will be estimated using the model defined in \code{data.rct$psName}.
#' @param rct.trControl An object of [caret::trainControl()] for [caret::train()], 
#'   which controls the computational nuances for fitting the outcome model 
#'   using the RCT dataset.
#' @param ec.trControl An object of [caret::trainControl()], which controls the
#'   computational nuances for fitting the outcome model using each EC dataset.
#' @param method A character object. The classification or regression model 
#'   [caret::train()] will use to estimate the outcome model.
#' @param nu.vector NULL or a numeric vector. The proposed nu values, the selected
#'   value minimizes the trade-off between...
#' @param min.lambda NULL, scalar numeric, or numeric vector.  The minimum value 
#'   of lambda to consider in the glmnet regression.  If NULL, 
#'   \eqn{min.lambda = n^-(1.0 + nu / 5.0) / 2.0}, where n is the total number
#'   of participants across the EC datasets and nu is each value of nu.vector.
#'   If a scalar, the same value is used for all values of nu.vector. If
#'   a numeric vector, it must be of the same length as nu.vector. The max is set as 
#'   \eqn{min(n^(-0.1), 0.9)}; 
#'   100 evenly spaced lambda values between these extrema are considered.
#'
#' @returns A list with components:
#' * est: estimated average treatment effect by AIPW, ACW and the selective integrative estimator.
#' * sd: estimated standard errors for the aforementioned estimators.
#' * subset.idx: a subset of indices of the external controls which have been selected for the
#' final integrative estimation.
#' 
#' @examples
#' data(selectiveToy)
#' 
#' data_rct <- dataInput(selectiveToy.rct, Y~A*X1 + A*X2, A~X1+X2)
#' 
#' # can manually construct data_ec
#' data_ec <- list("X" = data.matrix(selectiveToy.rwe[, c("X1","X2")]),
#'                 "Y" = selectiveToy.rwe$Y,
#'                 "A" = selectiveToy.rwe$A)
#'                 
#' # or use dataInput() to construct data.ec input using the rct models to
#' # ensure that all required covariates are kept.
#' data_ec <- dataInput(selectiveToy.rwe, Y~A*X1 + A*X2, A ~ X1+X2)
#' 
#' result <- srEC(data.rct = data_rct,
#'                data.ec = data_ec,
#'                method = "glm")
#' 
#' @include stopTests.R fitModels.R estimateAIPW.R
#' @export
srEC <- function(data.rct,
                 data.ec = NULL,
                 ...,
                 ps.rct = NULL,
                 rct.trControl = caret::trainControl(method = 'cv', number = 10L),
                 ec.trControl = caret::trainControl(method = 'cv', number = 10L),
                 method = "gbm",
                 nu.vector = c(1.0, 2.0),
                 min.lambda = NULL) {
  
  stopifnot(
    "`data.rct` must be provided" = !missing(data.rct)
  )
  
  # ensure that provided data match expected structure
  .isDI(data.rct, "data.rct")

  stopifnot(
    "`data.ec` must be a list" = is.null(data.ec) || is.vector(data.ec, "list"),
    "`ps.rct` must be NULL or a numeric vector of length = nrow(data.rct$X)" =
      is.null(ps.rct) || .isNumericVector(ps.rct, nrow(data.rct$X)),
    "`rct.trControl` must be a list()" = is.vector(rct.trControl, "list"),
    "`ec.trControl` must be a list()" = is.vector(ec.trControl, "list"),
    "`method` must be a character" = .isCharacterVector(method, 1L),
    "`nu.vector` must be a numeric vector" = .isNumericVector(nu.vector),
    "`min.lambda` must be NULL, scalar numeric or a numeric vector" =
      is.null(min.lambda) || .isNumericVector(min.lambda, 1L) ||
      .isNumericVector(min.lambda, length(nu.vector))
  )  
  
  n_rct <- length(data.rct$Y)
  
  data.rct$mainName <- .adjustModelCoding(data.rct$mainName, colnames(data.rct$X))
  data.rct$contName <- .adjustModelCoding(data.rct$contName, colnames(data.rct$X))
  data.rct$psName <- .adjustModelCoding(data.rct$psName, colnames(data.rct$X))
  
  if (!is.null(data.ec)) {
    if (all(c("X", "Y") %in% names(data.ec))) data.ec <- list(data.ec)
  }

  # fit the propensity score and outcome models using the RCT and EC
  # datasets, updating the respective lists to contain the fitted values
  # note that this procedure turns each element of data.ec into a list
  # equivalent to that generated by dataInputs(), augmented with the
  # fitted values
  updated_data_objects <- .fitModels(data.rct = data.rct, data.ec = data.ec,
                                     ps.rct = ps.rct,
                                     rct.trControl = rct.trControl,
                                     ec.trControl = ec.trControl,
                                     method = method, ...)
  rm(data.rct, data.ec)
  
  aipw_result <- .estimateAIPW(data.rct = updated_data_objects$data.rct)
  
  if (is.null(updated_data_objects$data.ec) || 
      length(updated_data_objects$data.ec) == 0L) {
    return( list("AIPW" = aipw_result[c("tau.hat", "sd.hat", "CI")]) )
  }

  acw_result <- .estimateACW(data.rct = updated_data_objects$data.rct, 
                             data.ec = updated_data_objects$data.ec,
                             aipw.result = aipw_result)
  
  bias_lasso <- .estimateLASSO(acw.result = acw_result, n.rct = n_rct,
                               nu.vector = nu.vector, min.lambda = min.lambda)
    
  # obtain the indices for borrowing
  ec_idx_lasso <- {abs(unlist(bias_lasso)) < 1e-8} |> which() |> unname()

  if (length(ec_idx_lasso) == 0L) {
    # if no external controls are selected
    return( list("aipw" = aipw_result[c("tau.hat", "sd.hat", "CI")],
                 "acw" = list("tau.hat" = acw_result$tau.hat,
                              "sd.hat" = acw_result$sd.hat),
                 "subset.idx" = integer(0L)) )
  }
  acw_lasso_result <- .estimateACWLASSO(data.rct = updated_data_objects$data.rct,
                                        data.ec = updated_data_objects$data.ec,
                                        aipw.result = aipw_result,
                                        bias.lasso = bias_lasso)
  
  acw_final <- .estimateFinal(aipw.result = aipw_result, 
                              acw.lasso.result = acw_lasso_result)
  
  list("aipw" = aipw_result[c("tau.hat", "sd.hat")],
       "acw" = acw_result[c("tau.hat", "sd.hat")],
       "acw.lasso" = acw_lasso_result[c("tau.hat", "sd.hat")],
       "acw.final" = acw_final,
       "subset.idx" = ec_idx_lasso)
}
