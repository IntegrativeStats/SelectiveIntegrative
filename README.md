# SelectiveIntegrative
Selective Integrative Analyses of the Average Treatment Effect Combining External Controls with the Standard Randomized Clinic Trial.

❗❗❗This is a late stage beta release. If you encounter any issues or have suggested changes, please reach out to shannon.holloway at duke.edu.


<h4>Usage</h4>
<pre>
srEC(data.rct, data.ec = NULL,
     ...,
     ps.rct = NULL,
     rct.trControl = caret::trainControl(method = "cv", number = 10L),
     ec.trControl = caret::trainControl(method = "cv", number = 10L),
     method = "gbm",
     nu.vector = c(1, 2),
     min.lambda = NULL
)
</pre>

<h5>Abbreviated Formal Argument Descriptions</h5>

- data.rct: The value object returned by *dataInput()* for the
  data from a randomized clinical trial (RCT). 
- data.ec: The value object returned by *dataInput()* for the
  data from an external control (EC). If not provided, only AIPW
  estimates will be returned.
- ...: Other options used to fit the predictive models. Passed on to
  caret::train().
- ps.rct: Optional input providing a vector of known propensity
  scores P(A=1) for the RCT dataset.
- rct.trControl: The value object returned by caret::trainControl()} to 
  controls the computational nuances for fitting the outcome model
  using the RCT dataset.
- ec.trControl: The value object returned by caret::trainControl()} to 
  controls the computational nuances for fitting the outcome model
  using the EC dataset.
- method: The classification or regression model caret::train() will use to 
  estimate the outcome model.
- nu.vector: The proposed nu values.
- min.lambda: The minimum value of lambda to consider in the glmnet regression.

<h5>Abbreviate Returned Object Description</h5>

A list with components:

- est: estimated average treatment effect by AIPW, ACW and the selective integrative estimator.
- sd: estimated standard errors for the aforementioned estimators.
- subset.idx: a subset of indices of the external controls which have been selected for the
final integrative estimation.

<h3>Examples</h3>

<pre>
  
# load provided illustrative toy dataset
data(selectiveToy)

data_rct <- dataInput(selectiveToy.rct, Y~A*X1 + A*X2, A~X1+X2)

# can manually construct data_ec
data_ec <- list("X" = data.matrix(selectiveToy.rwe[, c("X1","X2")]),
                "Y" = selectiveToy.rwe$Y,
                "A" = selectiveToy.rwe$A)
                
# or use dataInput() to construct data.ec input using the rct models to
# ensure that all required covariates are kept.
data_ec <- dataInput(selectiveToy.rwe, Y~A*X1 + A*X2, A ~ X1+X2)

result <- srEC(data.rct = data_rct,
               data.ec = data_ec,
               method = "glm")

# Expected Results:
cat("AIPW:", format(result$aipw$tau.hat, digits = 3), 
    "S.E.: ", format(result$aipw$sd.hat / sqrt(nrow(selectiveToy.rct)), digits = 3))
## AIPW: -0.229 S.E.:  0.161

cat("ACW:", format(result$acw$tau.hat, digits = 3), 
    "S.E.: ", format(result$acw$sd.hat / sqrt(nrow(selectiveToy.rct)), digits = 3))
## ACW: -0.706 S.E.:  0.162

cat("ACW.final:", format(result$acw.final$tau.hat, digits = 3), 
    "S.E.: ", format(result$acw.final$sd.hat / sqrt(nrow(selectiveToy.rct)), digits = 3))
## ACW.final: -0.235 S.E.:  0.159

</pre>

<h3>Install</h3>

In your R console

<pre>
library(remotes)
install_github("IntegrativeStats/SelectiveIntegrative")
</pre>
