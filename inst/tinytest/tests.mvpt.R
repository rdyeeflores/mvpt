#### COMPONENTS: Libraries, data, and EXAMPLEs for testing ####

library(lavaan) 
## library(tinytest)
## library(mvpt)

## DATA: Simulated dataset (no regression needed)
set.seed(12345)
data_generating_input <- 
  "
  ## latent variables
  LV1 =~ 0.5*x1  + 0.7*x2  + 0.9*x3
  LV2 =~ 5.0*x4  + 7.0*x5  + 9.0*x6
  LV3 =~ 0.5*x7  + 0.7*x8  + 0.9*x9
  "
simData <- lavaan::simulateData(model = data_generating_input,  sample.nobs = 300)
## Adding one NA for missingness
simData[1, 1] <- NA
## Adding composite measures
simData$C1 <- (simData$x1  + simData$x2  + simData$x3)  / 3
simData$C2 <- (simData$x4  + simData$x5  + simData$x6)  / 3
simData$C3 <- (simData$x7  + simData$x8  + simData$x9)  / 3

## EXAMPLE: One correct model wo LVs
corr_input <- 
  "
  C2 ~ C1
  C3 ~ C1 + C2
  "

## EXAMPLE: Sextet of LV models (comprise one whole MEC)
model1 <- 
  "
  ## regressions
  LV2 ~ LV1
  LV3 ~ LV1 + LV2
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  "
model2 <- 
  "
  ## regressions
  LV1 ~ LV2
  LV3 ~ LV1 + LV2
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  "
model3 <- 
  "
  ## regressions
  LV3 ~ LV1
  LV2 ~ LV1 + LV3
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  "
model4 <- 
  "
  ## regressions
  LV1 ~ LV3
  LV2 ~ LV1 + LV3
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  "
model5 <- 
  "
  ## regressions
  LV3 ~ LV2
  LV1 ~ LV2 + LV3
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  "
model6 <- 
  "
  ## regressions
  LV2 ~ LV3
  LV1 ~ LV2 + LV3
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  "

## EXAMPLE: Sextet of fitted LV models
model1_fit <- sem(model1, data=simData, fixed.x = FALSE, missing = "ml")
model2_fit <- sem(model2, data=simData, fixed.x = FALSE, missing = "ml")
model3_fit <- sem(model3, data=simData, fixed.x = FALSE, missing = "ml")
model4_fit <- sem(model4, data=simData, fixed.x = FALSE, missing = "ml")
model5_fit <- sem(model5, data=simData, fixed.x = FALSE, missing = "ml")
model6_fit <- sem(model6, data=simData, fixed.x = FALSE, missing = "ml")

## EXAMPLE: mvpt() output
mvpt_output <- mvpt(lavaan_input = model1, path = "LV3 ~ LV1", data = simData, showplots = TRUE)


#### THEORY TESTS (maybe add one more about automatics covs) ####

## Functionality of equal fit indices theory across SEXTET (within tolerance)
model_fit_list <- list(model1_fit, model2_fit, model3_fit, model4_fit, model5_fit, model6_fit)
fm_list <- lapply(model_fit_list, function(f) unname(fitMeasures(f, c("cfi","tli","rmsea","srmr"))))
expect_true(all(sapply(fm_list, function(x) isTRUE(all.equal(fm_list[[1]], x, tolerance = 1e-6)))))


#### MAIN FUNCTION: mvpt(lavaan_input, path, data) ####

## Functionality with LV model
expect_silent(mvpt_output)

## Functionality with LV model, list()
expect_silent(mvpt(lavaan_input = list(model1, model2), path = "LV3 ~ LV1", data = simData, showplots = TRUE))

## Functionality without LVs
expect_silent(mvpt(lavaan_input = corr_input, path = "C3 ~ C1", data = simData, showplots = TRUE))

## Functionality without LVs, list() with one model nested in the next (where path match is implied)
expect_silent(mvpt(lavaan_input = list("C3 ~ C1", "C3 ~ C2 + C1\n C2~C1"), path = "C3 ~ C1", data = simData, showplots = TRUE))

## Non-functionality due to negative variances (Heywood cases)
set.seed(123); n <- 100
f1 <- rnorm(n)
f2 <- 0.80 * f1 + rnorm(n, sd = 0.30)
y1 <- 1.00 * f1 + rnorm(n, sd = 0.01) ## perfect indicator
y2 <- 0.40 * f1 + rnorm(n, sd = 0.30)
y3 <- 0.40 * f2 + rnorm(n, sd = 0.30)
y4 <- 0.40 * f2 + rnorm(n, sd = 0.30)
dat_hey <- data.frame(y1, y2, y3, y4)
model_hey1 <- 
  '
  f1 =~ y1 + y2
  f2 =~ y3 + y4
  f2 ~ f1
  '
model_hey2 <- 
  '
  f1 =~ y1 + y3
  f2 =~ y2 + y4
  f2 ~ f1
  '
lavaan_input <- list(model_hey1, model_hey2)
expect_error(mvpt(lavaan_input, path = "f2 ~ f1", data = dat_hey))

## Non-functionality due to non-positive definite cov matrix
set.seed(123); n <- 100
x1 <- rnorm(n)
x2 <- x1          # perfect correlation
x3 <- rnorm(n)
y  <- rnorm(n)
data <- data.frame(x1, x2, x3, y)
lavaan_input <- 
  '
  y ~ x1 + x2
  x2 ~ x1
  '
path = "y ~ x1"
expect_error(mvpt(lavaan_input, path = "y ~ x1", data))

## Non-functionality due wrong input object
expect_error(mvpt(lavaan_input = model1_fit, path = "C3~C1", data = simData))

## Non-functionality due non-lavaan formatting error
err_input <- 
  "
  C2 + C1
  C3 ~ C1 + C2
  "
expect_error(mvpt(lavaan_input = list(corr_input, err_input), path = "C3~C1", data = simData, showplots = TRUE))

## Non-functionality message due to covariance specification (update pending) 
err_input <- 
  "
  C2 ~~ C1
  C3 ~ C1 + C2
  "
expect_error(mvpt(lavaan_input = list(corr_input, err_input), path = "C3~C1", data = simData, showplots = TRUE))

## Non-functionality due to label usage
err_input <- 
  "
  C2 ~ a*C1
  C3 ~ C1 + C2
  "
expect_error(mvpt(lavaan_input = list(corr_input, err_input), path = "C3~C1", data = simData, showplots = TRUE))

## Non-functionality due to path misspecification
expect_error(mvpt(lavaan_input = corr_input, path = "C3->C1", data = simData))

## Non-functionality due to path not being in model
expect_error(mvpt(lavaan_input = corr_input, path = "C1~C3", data = simData))

## Non-functionality due to list() input and specifying a path that is not shared
expect_error(mvpt(lavaan_input = list(model1, model4), path = "LV3 ~ LV1", data = simData, showplots = TRUE))


#### MAIN FUNCTION: mvptRank(mvpt_output) ####

## Functionality assuming mvpt() functionality
expect_silent(mvptRank(mvpt_output))

## Non-functionality assuming non mvpt class
expect_error(mvptRank(model1_fit))


#### MAIN FUNCTION: mvptZoom(mvpt_output, index) ####

## Functionality assuming mvpt() functionality
expect_silent(mvptZoom(mvpt_output, index=1))

## Non-functionality assuming non mvpt class
expect_error(mvptZoom(model1_fit, index=1))

## Non-functionality due to out of bounds index
expect_error(mvptZoom(mvpt_output, index=7))


#### HELPER FUNCTION: dagu(LAV_list, path) ####

LAV_list <- "LV3 ~ LV1"
path <- "LV3 ~ LV1"

## Non-functionality message if user-given model is orphan
expect_message(mvpt:::dagu(LAV_list, path))


#### HELPER FUNCTION: auto_sem(fam_lavaan_ready, data, missing = "ml", estimator = "ML") ####

fam_lavaan_ready <- mvpt:::dagu(LAV_list = model1, path = "LV3 ~ LV1")$fam_lavaan_ready
data <- simData

## Non-functionality if default missing argument is changed
expect_error(auto_sem(fam_lavaan_ready, data, missing = "listwise"))

## Non-functionality if default estimator argument is changed
expect_error(auto_sem(fam_lavaan_ready, data, estimator = "WLS"))


#### HELPER FUNCTION: calc_A(SEMfitted), calc_sc(SEMfitted), & calc_B(sc1, sc2, n) ####

SEMfitted <- model1_fit
sc1 <- mvpt:::calc_sc(SEMfitted)
sc2 <- mvpt:::calc_sc(SEMfitted)   
n <- nobs(SEMfitted)


#### HELPER FUNCTION: vw_core(SEMfitted_list, path) ####

SEMfitted_list <- list(model1_fit, model2_fit, model3_fit, model4_fit, model5_fit)
path <- "LV3 ~ LV1"





