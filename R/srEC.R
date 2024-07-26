

#' Selective and robust external-control borrowing strategy for evaluating treatment effect
#'
#' @description
#' `srEC()` is a penalized dynamic integrative framework to augment a randomized trial (RT)
#'  with a external control (EC) dataset, in which the subject-level compatibility of the EC
#'  is assessed by a well-crafted penalty (e.g., adaptive lasso penalty). The parameter of interest
#'  is the average treatment effect.
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom caret train
#' @param data.rct A list contains X, Y and A for RT. The propensity scores P(A=1|X) can also be
#' contained as `ps` (or `NULL`).
#' @param data.ec A list contains X and Y for EC with treatment A=0.
#' @param rct.trControl An object of [caret::trainControl()] for [caret::train()], which controls the
#' computational nuances for fitting the outcome model using the RT dataset.
#' @param ec.trControl An object of [caret::trainControl()], which controls the
#' computational nuances for fitting the outcome model using the EC dataset.
#' @param ... Other options used to fit the predictive models. Passed on to [caret::train()].
#'
#' @returns A list with components:
#' * est: estimated average treatment effect by AIPW, ACW and the selective integrative estimator.
#' * sd: estimated standard errors for the aforementioned estimators.
#' * subset.idx: a subset of indices of the external controls which have been selected for the
#' final integrative estimation.
#' @export
srEC <- function(data.rct,
                 data.ec = NULL,
                 rct.trControl = caret::trainControl(method = 'cv', number = 10),
                 ec.trControl = caret::trainControl(method = 'cv', number = 10),
                 method = 'gbm', ...) {
  
  # assign the variables for analyses
  n_rct <- length(data.rct$Y)

  updated_data_objects <- .fitModels(data.rct = data.rct, data.ec = data.ec,
                                     rct.trControl = rct.trControl,
                                     ec.trControl = ec.trControl)
  rm(data.rct, data.ec)
  
  aipw_result <- .aipw(data.rct = updated_data_objects$data.rct)
  
  if (is.null(updated_data.ec) || length(updated_data.ec) == 0L) {
    return(list("AIPW" = apiw_result[c("tau.hat", "sd.hat", "CI")]))
  }

  acw_result <- .estimateACW(data.rct = updated_data_objects$data.rct, 
                             data.ec = updated_data_objects$data.ec,
                             aipw.result = aipw_result)
  
  bias_lasso <- .estimateLASSO(acw.result = acw.result, 
                               n.rct = n_rct)
    
  # obtain the indices for borrowing
  ec_idx_lasso <- unname(which(abs(unlist(bias_lasso)) < 1e-8))

  if (length(ec_idx_lasso) == 0L) {
    # if no external controls are selected
    return(list("aipw" = aipw_result[c("tau.hat", "sd.hat", "CI")],
                "acw" = list("tau.hat" = acw_result$tau.hat,
                             "sd.hat" = acw_result$sd.hat),
                "subset.idx" = integer(0L)))
  }
  
  acw_lasso_result <- .estimateACWLASSO(data.rct = data_objects$data.rct,
                                        data.ec = data_objects$data.ec,
                                        bias.lasso = bias_lasso)
  
  acw_final <- .estimateFinal(aipw.result = aipw_result, 
                              acw.lasso.result = acw_lasso_result)
  
  return(list("aipw" = aipw_result[c("tau.hat", "sd.hat")],
              "acw" = acw_result[c("tau.hat", "sd.hat")],
              "acw.lasso" = acw_lasso_result[c("tau.hat", "sd.hat")],
              "acw.final" = acw_final,
              "subset.idx" = ec_idx_lasso))
}
