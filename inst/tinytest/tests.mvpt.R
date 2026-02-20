#### COMPONENTS: Libraries, data, and fitted models ####

library(lavaan)
library(dagitty)

## DATA: Simulated workburnout dataset (no regression yet)
set.seed(12345)
data_generating_model <- 
  "
  ## latent variables
  LV1 =~ 0.5*x1  + 0.7*x2  + 0.9*x3
  LV2 =~ 0.5*x4  + 0.7*x5  + 0.9*x6
  LV3 =~ 0.5*x7  + 0.7*x8  + 0.9*x9
  LV4 =~ 0.5*x10 + 0.7*x11 + 0.9*x12
  LV5 =~ 0.5*x13 + 0.7*x14 + 0.9*x15
  LV6 =~ 0.5*x16 + 0.7*x17 + 0.9*x18
  ## regressions (none)
  "
simData <- simulateData(model = data_generating_model,  sample.nobs = 599)
## Adding one NA to this complete dataset
simData[10, 10] <- NA
## Adding composite measures
simData$C1 <- (simData$x1  + simData$x2  + simData$x3)  / 3
simData$C2 <- (simData$x4  + simData$x5  + simData$x6)  / 3
simData$C3 <- (simData$x7  + simData$x8  + simData$x9)  / 3
simData$C4 <- (simData$x10 + simData$x11 + simData$x12) / 3
simData$C5 <- (simData$x13 + simData$x14 + simData$x15) / 3
simData$C6 <- (simData$x16 + simData$x17 + simData$x18) / 3
## Five LV models (all members of same MEC)
simM1 <- 
  "
  ## regressions
  LV4 ~ LV1
  LV1 ~ LV6
  LV4 ~ LV3
  LV3 ~ LV2
  LV6 ~ LV5
  LV4 ~ LV5
  LV2 ~ LV1
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  LV4 =~ x10 + x11 + x12
  LV5 =~ x13 + x14 + x15
  LV6 =~ x16 + x17 + x18
  "
simM2 <- 
  "
  ## regressions
  LV4 ~ LV1
  LV5 ~ LV6
  LV1 ~ LV6
  LV4 ~ LV3
  LV3 ~ LV2
  LV4 ~ LV5
  LV2 ~ LV1
## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  LV4 =~ x10 + x11 + x12
  LV5 =~ x13 + x14 + x15
  LV6 =~ x16 + x17 + x18
  "
simM3 <- 
  "
  ## regressions
  LV4 ~ LV1
  LV5 ~ LV6
  LV4 ~ LV3
  LV3 ~ LV2
  LV4 ~ LV5
  LV6 ~ LV1
  LV2 ~ LV1
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  LV4 =~ x10 + x11 + x12
  LV5 =~ x13 + x14 + x15
  LV6 =~ x16 + x17 + x18
  "
simM4 <- 
  "
  ## regressions
  LV4 ~ LV1
  LV5 ~ LV6
  LV4 ~ LV3
  LV3 ~ LV2
  LV1 ~ LV2
  LV4 ~ LV5
  LV6 ~ LV1
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  LV4 =~ x10 + x11 + x12
  LV5 =~ x13 + x14 + x15
  LV6 =~ x16 + x17 + x18
  "
simM5 <- 
  "
  ## regressions
  LV4 ~ LV1
  LV5 ~ LV6
  LV2 ~ LV3
  LV4 ~ LV3
  LV1 ~ LV2
  LV4 ~ LV5
  LV6 ~ LV1
  ## latent variables
  LV1 =~ x1  + x2  + x3
  LV2 =~ x4  + x5  + x6
  LV3 =~ x7  + x8  + x9
  LV4 =~ x10 + x11 + x12
  LV5 =~ x13 + x14 + x15
  LV6 =~ x16 + x17 + x18
  "
simM1_fit <- sem(simM1, data=simData, fixed.x = FALSE, missing = "ml")
simM2_fit <- sem(simM2, data=simData, fixed.x = FALSE, missing = "ml")
simM3_fit <- sem(simM3, data=simData, fixed.x = FALSE, missing = "ml")
simM4_fit <- sem(simM4, data=simData, fixed.x = FALSE, missing = "ml")
simM5_fit <- sem(simM5, data=simData, fixed.x = FALSE, missing = "ml")


#### MAIN FUNCTION: mvpt(lavaan_model, path, data) ####

## Functionality with LVs
mvpt_output <- mvpt(lavaan_model = simM1, path = "LV4 ~ LV1", data = simData, showplots = TRUE)
expect_silent(mvpt_output)

## Functionality without LVs
corr_model_format <- 
"
  C4 ~ C1
  C5 ~ C6
  C4 ~ C3
  C3 ~ C2
  C1 ~ C2
  C4 ~ C5
  C6 ~ C1
"
expect_silent(mvpt(lavaan_model = corr_model_format, path = "C4 ~ C1", data = simData, showplots = TRUE))

## Non-functionality due to label usage
err_model_format <- 
"
  C4 ~ a*C1
  C5 ~ C6
  C4 ~ C3
  C3 ~ C2
  C1 ~ C2
  C4 ~ C5
  C6 ~ C1
"
expect_error(mvpt(lavaan_model = err_model_format, path = "C4~C1", data = simData, showplots = TRUE))

## Non-functionality due non-lavaan formatting 
err_model_format <- 
"
  C4 + C1
  C5 ~ C6
  C4 ~ C3
  C3 ~ C2
  C1 ~ C2
  C4 ~ C5
  C6 ~ C1
"
expect_error(mvpt(lavaan_model = err_model_format, path = "C4~C1", data = simData, showplots = TRUE))

## Non-functionality message due to covariance specification (update pending) 
err_model_format <- 
  "
  C4 ~~ C1
  C5 ~ C6
  C4 ~ C3
  C3 ~ C2
  C1 ~ C2
  C4 ~ C5
  C6 ~ C1
"
expect_error(mvpt(lavaan_model = err_model_format, path = "C4~C1", data = simData, showplots = TRUE))

## Non-functionality and warnings due to lavaan fitting error (!need a list first)


#### MAIN FUNCTION: mvptRank(mvpt_output) ####

## Functionality assuming mvpt() functionality
expect_silent(mvptRank(mvpt_output))


#### MAIN FUNCTION: mvptZoom(mvpt_output, index) ####

## Functionality assuming mvpt() functionality
expect_silent(mvptZoom(mvpt_output, index=1))

## Non-functionality due to out of bounds index
expect_error(mvptZoom(mvpt_output, index=7))


#### HELPER FUNCTIONS ####

LAV <- simM1
path <- "LV4 ~ LV1"
data <- simData
subMEC_lavaan_ready <- mvpt:::dagu(LAV, path)$subMEC_lavaan_ready
mvpt:::auto_sem(subMEC_lavaan_ready, data)

## Non-functionality if default missing argument is changed
expect_error(auto_sem(subMEC_lavaan_ready, data, missing = "listwise"))

## Non-functionality if default estimator argument is changed
expect_error(auto_sem(subMEC_lavaan_ready, data, estimator = "WLS"))

SEMfitted <- simM1_fit
sc1 <- mvpt:::calc_sc(SEMfitted)
sc2 <- mvpt:::calc_sc(SEMfitted)   
n <- nobs(SEMfitted)
mvpt:::calc_A(SEMfitted)
mvpt:::calc_sc(SEMfitted)
mvpt:::calc_B(sc1, sc2, n)

SEMfitted_list <- list(simM1_fit, simM2_fit, simM3_fit, simM4_fit, simM5_fit)
path <- "LV4~LV1"
mvpt:::VW_core(SEMfitted_list, path)

