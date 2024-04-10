## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE
)

## ----setup--------------------------------------------------------------------
# devtools::install_github("yangzhao98/drugTargetScreen",force=TRUE)
library(drugTargetScreen)
path.res <- getwd()

## ----eval=FALSE---------------------------------------------------------------
#  ## preparation
#  source(paste(path.res,"/0_precondition.R",sep=""))
#  source(paste(path.res,"/utilities.R",sep=""))
#  
#  ## codes for stage 1. Causal biomarkers screen
#  source(paste(path.res,"/1_cause_biomarkers_screen.R",sep=""))
#  
#  ## codes for stage 2. Drug-targeting protein discovery and verification
#  source(paste(path.res,"/2_discovery_verification.R",sep=""))
#  
#  ## codes for stage 3. Drug target relevant side effects exploration
#  source(paste(path.res,"/3_side_effects.R",sep=""))
#  
#  ## codes for Stage 4. Therapeutic efficacy and safety emulation
#  source(paste(path.res,"/4_target_trial_emulation.R",sep=""))

## -----------------------------------------------------------------------------
load(paste(
  path.res,"/",
  grep(".RData",list.files(path.res),value=TRUE),
  sep=""))

## ----fig.width=7,fig.height=5-------------------------------------------------
plt.woCor

## ----fig.width=7,fig.height=5-------------------------------------------------
plt.4Gene

## ----fig.width=7,fig.height=5-------------------------------------------------
plt.GRAPPLE.profile

## ----fig.width=7,fig.height=5-------------------------------------------------
plt.Forest

## ----fig.width=20,fig.height=10-----------------------------------------------
plt.Step2

## ----fig.width=20,fig.height=5------------------------------------------------
plt.Step3

## ----fig.width=20,fig.height=5------------------------------------------------
plt.Step4

