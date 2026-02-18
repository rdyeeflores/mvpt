
#### COMPONENTS: Libraries, data, and fitted models ####

## pkgload::load_all(export_all = FALSE) ## warning expected

library(mvpt)
library(lavaan)
library(dagitty)

## DATA: Workburnout 
load("C:/Git/mvpt-research/DATA/burnout/burnout.rda")
## Add one NA to this complete dataset
burnout[200, 15] <- NA
anyNA(burnout)
## Adding composites
burnout$DMc <-  (burnout$DM1 + burnout$DM2) / 2
burnout$SEc <-  (burnout$SE1 + burnout$SE2 + burnout$SE3) / 3
burnout$ELCc <- (burnout$ELC1 + burnout$ELC2 + burnout$ELC3 + burnout$ELC4 + burnout$ELC5) / 5
burnout$EEc <-  (burnout$EE1 + burnout$EE2 + burnout$EE3) / 3
burnout$DPc <-  (burnout$DP1 + burnout$DP2) / 2
burnout$PAc <-  (burnout$PA1 + burnout$PA2 + burnout$PA3) / 3

## Five LV models (all members of same MEC) 
M1 <- 
  "
  ## regressions
  SE ~ DM
  PA ~ DP
  DP ~ EE
  DM ~ ELC
  PA ~ ELC
  EE ~ SE
  PA ~ SE
  ## latent variables
  DM  =~ DM1 + DM2
  SE  =~ SE1 + SE2 + SE3
  ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
  EE  =~ EE1 + EE2 + EE3
  DP  =~ DP1 + DP2
  PA  =~ PA1 + PA2 + PA3
  "
M2 <- 
  "
  ## regressions
  ELC ~ DM
  SE ~ DM
  PA ~ DP
  DP ~ EE
  PA ~ ELC
  EE ~ SE
  PA ~ SE
## latent variables
  DM  =~ DM1 + DM2
  SE  =~ SE1 + SE2 + SE3
  ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
  EE  =~ EE1 + EE2 + EE3
  DP  =~ DP1 + DP2
  PA  =~ PA1 + PA2 + PA3
  "
M3 <- 
  "
  ## regressions
  ELC ~ DM
  PA ~ DP
  DP ~ EE
  PA ~ ELC
  DM ~ SE
  EE ~ SE
  PA ~ SE
  ## latent variables
  DM  =~ DM1 + DM2
  SE  =~ SE1 + SE2 + SE3
  ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
  EE  =~ EE1 + EE2 + EE3
  DP  =~ DP1 + DP2
  PA  =~ PA1 + PA2 + PA3
  "
M4 <- 
  "
  ## regressions
  ELC ~ DM
  PA ~ DP
  DP ~ EE
  SE ~ EE
  PA ~ ELC
  DM ~ SE
  PA ~ SE
  ## latent variables
  DM  =~ DM1 + DM2
  SE  =~ SE1 + SE2 + SE3
  ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
  EE  =~ EE1 + EE2 + EE3
  DP  =~ DP1 + DP2
  PA  =~ PA1 + PA2 + PA3
  "
M5 <- 
  "
  ## regressions
  ELC ~ DM
  EE ~ DP
  PA ~ DP
  SE ~ EE
  PA ~ ELC
  DM ~ SE
  PA ~ SE
  ## latent variables
  DM  =~ DM1 + DM2
  SE  =~ SE1 + SE2 + SE3
  ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
  EE  =~ EE1 + EE2 + EE3
  DP  =~ DP1 + DP2
  PA  =~ PA1 + PA2 + PA3
  "
M1_fit <- sem(M1, data=burnout, fixed.x = FALSE)
M2_fit <- sem(M2, data=burnout, fixed.x = FALSE)
M3_fit <- sem(M3, data=burnout, fixed.x = FALSE)
M4_fit <- sem(M4, data=burnout, fixed.x = FALSE)
M5_fit <- sem(M5, data=burnout, fixed.x = FALSE)


#### FUNCTION: dagu(LAV, path) and auto_sem(subMEC_lavaan_ready, data) ####

LAV <-  ## correct wo LVs
  "
  PAc ~ EEc + DPc
  DPc ~ EEc
  "
LAV <-   ## correct w LVs
  "
  ## regs
  PA ~ EE + DP
  DP ~ EE
  ## LVs
  EE  =~ EE1 + EE2 + EE3
  DP  =~ DP1 + DP2
  PA  =~ PA1 + PA2 + PA3
  "
LAV <- ## incorrect: covariance present
  "
  PAc ~ EEc + DPc
  DPc ~~ EEc
  "
LAV <- ## incorrect: adding parameter labels
  "
  PAc ~ a*EEc + c*DPc
  DPc ~ b*EEc
  "
path <- "PAc~EEc" ## correct wo LVs
path <- "PA~EE"   ## correct w LVs
path <- "PA ~ EE" ## correct w LVs but with gaps around ~
dagu(LAV, path)
subMEC_lavaan_ready <- dagu(LAV, path)$subMEC_lavaan_ready
auto_sem(subMEC_lavaan_ready, data)
## Functionality without LVs
## Functionality with LVs
## Non-functionality with input formatting errors
## Non-functionality message with covariances 
## Non-functionality message with labels





#### FUNCTION: calc_A(SEMfitted), calc_sc(SEMfitted), and calc_B(sc1, sc2, n) #### 

SEMfitted <- M1_fit
sc1 <- calc_sc(M1_fit)
sc2 <- calc_sc(M2_fit)   
n <- nobs(M1_fit)
calc_A(SEMfitted)
calc_sc(SEMfitted)
calc_B(sc1, sc2, n)
## Functionality using fitted models
## Non-functionality because of library overlap (sandwich problem)
## Non-functionality because of matrix algebra error


#### FUNCTION: VW_core(SEMfitted_list, path) ####

SEMfitted_list <- list(M1_fit, M2_fit, M3_fit, M4_fit, M5_fit)
path <- "PA~SE"
path <- "PA ~ SE"
VW_core(SEMfitted_list, path)
## Functionality assuming all the above
## Functionality assuming all the above, but with lavaan warnings (eg: neg var)
## Non-functionality with impossible matrices




#### FUNCTION: mvpt(lavaan_model, path, data) and other main functions ####

lavaan_model <- 
  "
  ## regressions
  SE ~ DM
  PA ~ DP
  DP ~ EE
  DM ~ ELC
  PA ~ ELC
  EE ~ SE
  PA ~ SE
  ## latent variables
  DM  =~ DM1 + DM2
  SE  =~ SE1 + SE2 + SE3
  ELC =~ ELC1 + ELC2 + ELC3 + ELC4 + ELC5
  EE  =~ EE1 + EE2 + EE3
  DP  =~ DP1 + DP2
  PA  =~ PA1 + PA2 + PA3
  "
path <- "PA~SE"
path <- "PA ~ SE"
data <- burnout
MVP <- mvpt(lavaan_model, path, data, showplots = TRUE) ## 5 models; p=0.999
MVP
mvptZoom(MVP, 1)
## Functionality assuming all the above
## Non-functionality with input format errors





