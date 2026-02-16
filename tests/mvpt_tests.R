
#### mvpt.R: QUALITY CONTROL TESTS ##########################################

## NOTE: Good to automate these tests so one command can be used to get a report about what failed.

## NOTE: Below are copies of code chunks in mvpt.R -> the goal is to have a series of tests, such as:
## 1. Checking for valid lavaan syntax, producing messages ow; valid syntax for all used R packages
## 2. Checking for valid mvpt syntax, producing messages ow
## 3. Making sure model generation is as expected, and constrain where needed (eg: too many models generated)
## 4. Having messages for all possible VW_core test problems


########################################################################################################
#### CODE CHECKING  ####
########################################################################################################

if(FALSE){
  
  ## pkgload::load_all("../mvpt", export_all = FALSE) ## warning expected
  
  library(mvpt)
  library(lavaan)
  library(dagitty)
  
  
  
  ## DATA: Workburnout 
  load("DATA//burnout/burnout.rda")
  ## Adding composites
  burnout$DMc <-  (burnout$DM1 + burnout$DM2) / 2
  burnout$SEc <-  (burnout$SE1 + burnout$SE2 + burnout$SE3) / 3
  burnout$ELCc <- (burnout$ELC1 + burnout$ELC2 + burnout$ELC3 + burnout$ELC4 + burnout$ELC5) / 5
  burnout$EEc <-  (burnout$EE1 + burnout$EE2 + burnout$EE3) / 3
  burnout$DPc <-  (burnout$DP1 + burnout$DP2) / 2
  burnout$PAc <-  (burnout$PA1 + burnout$PA2 + burnout$PA3) / 3
  ## Add one NA to this complete dataset
  burnout[200, 15] <- NA
  anyNA(burnout)
  
  
  ## EXAMPLE: Five LV models fitted to the same data (all members of same MEC) 
  ## NOTE: Absence of fixed.x=FALSE during fitting
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
  M1_fit <- sem(M1, data=burnout)
  M2_fit <- sem(M2, data=burnout)
  M3_fit <- sem(M3, data=burnout)
  M4_fit <- sem(M4, data=burnout)
  M5_fit <- sem(M5, data=burnout)
  
  ## dagu(LAV, path) and auto_sem(subMEC_lavaan_ready, data) 
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
  LAV <- ## wrong: adding covariance
    "
    PAc ~ EEc + DPc
    DPc ~~ EEc
    "
  LAV <- ## wrong: adding parameter labels
    "
    PAc ~ a*EEc + c*DPc
    DPc ~ b*EEc
    "
  path <- "PAc~EEc"
  path <- "PA~EE"
  
  mvpt(LAV, path, data=burnout, showplots = TRUE)
  ## Looking at key component of auto_sem()
  subMEC_lavaan_ready <- dagu(LAV, path)$subMEC_lavaan_ready
  
  ## calc_A(SEMfitted) and calc_sc(SEMfitted) 
  SEMfitted <- M1_fit
  
  ## calc_B(sc1, sc2, n) 
  sc1 <- calc_sc(M1_fit)
  sc2 <- calc_sc(M2_fit)   
  n <- nobs(M1_fit)
  
  ## VW_core(SEMfitted_list, path) 
  SEMfitted_list <- list(M1_fit, M2_fit, M3_fit, M4_fit, M5_fit)
  path <- "PA~SE"
  
  ## mvpt(LAV, path, data) and mvptZoom(MVP, M) 
  data <- burnout
  path <- "PA~SE"
  LAV <- 
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
  MVP <- mvpt(LAV=LAV, path, data, showplots = TRUE) ## 5 models; p=0.999
  MVP
  mvptZoom(MVP, 1)
  mvpt(LAV=LAV, path="PA~ELC", data=burnout) ## 5 models; p=1.000 
  mvpt(LAV=LAV, path="DM~ELC", data=burnout) ## err: orphan model
  mvpt(LAV=LAV, path="SE~DM", data=burnout) ## 2 models; p=0.993
  mvpt(LAV=LAV, path="EE~SE", data=burnout) ## 3 models; p=0.998
  mvpt(LAV=LAV, path="DP~EE", data=burnout) ## 4 models; p=1.000
  mvpt(LAV=LAV, path="PA~DP", data=burnout) ## 5 models; p=1.000
  
  
}





