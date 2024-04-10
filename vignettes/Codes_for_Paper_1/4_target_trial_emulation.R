## intervention
library(data.table)
library(gfoRmula)

## --- Define hypothetical interventions
## per 1mmol/L LDL-C reduction -> reduce 21% CVD risk
ldlc_int_1mmol <- function(newdf,pool,intvar,intvals,time_name,t) {
    newdf[,`:=`(LDLC=pmax(get(intvar)-1,0),lipidDrugStatus=1)]
    # newdf[,`:=`(LDLC=pmax(get(intvar)-1,0))]
}

ldlc_int_statin <- function(newdf,pool,intvar,intvals,time_name,t) {
    # newdf[,`:=`(prop=runif(.N,0.3,0.49))]
    newdf[,`:=`(LDLC=get(intvar)*(1-0.4),lipidDrugStatus=1)]
    # newdf[,`:=`(LDLC=pmin(get(intvar),get(intvar)*(1-0.5)))]
}
int.list <- list(
    list(c(ldlc_int_1mmol)),
    list(c(ldlc_int_statin))
)
int.desc <- c(
    "per 1 mmol/L reduction",
    "40% reduction (Statins,HMGCR)"
)
intvar.list <- as.list(rep("LDLC",length(int.list)))
highrisk.id <- unique(datPP$ID[datPP$t0==0 & datPP$LDLC>=1.8])

## --- outcomes of interest
## model for outcome of interest
Ymodel <- function(outcomeY_name) {
    as.formula(
        paste(outcomeY_name,"~",
              paste(c(covnames_tf,covnames_cumavg,covnames_lag),collapse="+",sep=""),"+", 
              paste("cumavg_",covnames_int,collapse="+",sep=""),"+", 
              "t0+I(t0^2)+I(t0^3)",sep="" 
        )
    )
}
## All-cause mortality
outcomeY0_name  <- "AnyEvent"
outcomeY0_model <- Ymodel(outcomeY_name=outcomeY0_name)
outcomeY0_model 
## ASCVD
outcomeY1_name  <- "ASCVD"
outcomeY1_model <- Ymodel(outcomeY_name=outcomeY1_name)
outcomeY1_model
## CVD
outcomeY2_name  <- "CVD"
outcomeY2_model <- Ymodel(outcomeY_name=outcomeY2_name)
outcomeY2_model
## --- competing risks
## ASCVD
compeventD1_name  <- "ASCVDcp"
compeventD1_model <- Ymodel(outcomeY_name=compeventD1_name)
compeventD1_model
## CVD
compeventD2_name  <- "CVDcp"
compeventD2_model <- Ymodel(outcomeY_name=compeventD2_name)
compeventD2_model

## ----------------------------------------------------------
## pooled logistic regression model development
format95CI <- function(x,lci,uci) {
    paste(sprintf("%.2f",x), " (",
          sprintf("%.2f",lci), " to ",
          sprintf("%.2f",uci), ")", sep="")
}
getRiskRatio <- function(gform.fit,int.lab) {
    dat <- gform.fit$result
    dat$Int <- factor(dat$`Interv.`,levels=0:3,labels=int.lab)
    dat <- dat[,.SD[.N],by=.(Int)]
    dat <- dat[,.(`NP Risk`=sprintf("%.3f",`NP Risk`),
                  `g-form risk`=sprintf("%.3f",`g-form risk`),
                  RiskRatio=format95CI(`Risk ratio`,`RR lower 95% CI`,`RR upper 95% CI`),
                  RiskDiff=format95CI(`Risk difference`,`RD lower 95% CI`,`RD upper 95% CI`))
               ,by=.(Int)]
}
ncores <- parallel::detectCores()
nsimul <- 1e4

int.lab <- c(
  "No interention",
  "per 1 mmol/L LDL-C reduction",
  "40% reduction (Statins, HMGCR)",
  "61% reduction (PCSK9 inhibitors)")

# -- For all-cause mortality
system.time(
    gform_AllCause0 <- gfoRmula::gformula_survival(
        obs_data=datPP, id=id, time_points=time_points, time_name=time_name,
        basecovs=covnames_tf,
        covnames=covnames_tv, covtypes=covnames_tv_types,
        histories=histories_tv, histvars=histvars_tv,
        covparams=covparms_tv,
        outcome_name=outcomeY0_name, ymodel=outcomeY0_model,
        intvars=NULL, interventions=NULL, int_descript=NULL,
        nsimul=nsimul, seed=1234, parallel=TRUE, ncores=ncores
    )
)
print("Model0 for All-cause mortality")
print(gform_AllCause0)

# -- For ASCVD
system.time(
    gform_ASCVD0 <- gfoRmula::gformula_survival(
        obs_data=datPP, id=id, time_points=time_points, time_name=time_name,
        basecovs=covnames_tf,
        covnames=covnames_tv, covtypes=covnames_tv_types,
        histories=histories_tv, histvars=histvars_tv,
        covparams=covparms_tv,
        outcome_name=outcomeY1_name, ymodel=outcomeY1_model,
        compevent_name=compeventD1_name, compevent_model=compeventD1_model,
        intvars=NULL, interventions=NULL, int_descript=NULL,
        nsimul=nsimul, seed=1234, parallel=TRUE, ncores=ncores
    )
)
print("Model0 for ASCVD")
print(gform_ASCVD0)
# coef(gform_ASCVD0)

# -- for CVD
system.time(
    gform_CVD0 <- gfoRmula::gformula_survival(
        obs_data=datPP, id=id, time_points=time_points, time_name=time_name,
        basecovs=covnames_tf,
        covnames=covnames_tv, covtypes=covnames_tv_types,
        histories=histories_tv, histvars=histvars_tv,
        covparams=covparms_tv,
        outcome_name=outcomeY2_name, ymodel=outcomeY2_model,
        compevent_name=compeventD2_name, compevent_model=compeventD2_model,
        intvars=NULL, interventions=NULL, int_descript=NULL,
        nsimul=nsimul, seed=1234, parallel=TRUE, ncores=ncores
    )
)
print("Model0 for CVD")
print(gform_CVD0)
# coef(gform_CVD0)
save(gform_AllCause0,gform_ASCVD0,gform_CVD0,
     file=paste(path.res,"/pgform_Model0_Results.RData",sep="")
)

## Target trial emulation
## ---
## All available participants with baseline LDLC>=2.8mmol/L
system.time(
    gform_AllCause <- gfoRmula::gformula_survival(
        obs_data=datPP[datPP$ID %in% highrisk.id,], 
        id=id, time_points=time_points, time_name=time_name,
        basecovs=covnames_tf,
        covnames=covnames_tv, covtypes=covnames_tv_types,
        histories=histories_tv, histvars=histvars_tv,
        covparams=covparms_tv,
        outcome_name=outcomeY0_name, ymodel=outcomeY0_model,
        # compevent_name=compeventD1_name, compevent_model=compeventD1_model,
        intvars=intvar.list, interventions=int.list, int_descript=int.desc,
        nsimul=nsimul, seed=1234, nsamples=nboot, parallel=TRUE, ncores=ncores
    )
)
print("*** Hypothetical interventions for all-cause mortality")
print(gform_AllCause)
getRiskRatio(gform_AllCause,int.lab=int.lab)

# coef(gform_AllCause)
# -- for ASCVD
system.time(
    gform_ASCVD <- gfoRmula::gformula_survival(
        obs_data=datPP[datPP$ID %in% highrisk.id,], 
        id=id, time_points=time_points, time_name=time_name,
        basecovs=covnames_tf,
        covnames=covnames_tv, covtypes=covnames_tv_types,
        histories=histories_tv, histvars=histvars_tv,
        covparams=covparms_tv,
        outcome_name=outcomeY1_name, ymodel=outcomeY1_model,
        compevent_name=compeventD1_name, compevent_model=compeventD1_model,
        intvars=intvar.list, interventions=int.list, int_descript=int.desc,
        nsimul=nsimul, seed=1234, nsamples=nboot, parallel=TRUE, ncores=ncores
    )
)
print("*** Hypothetical interventions for ASCVD")
print(gform_ASCVD)
# coef(gform_ASCVD)
getRiskRatio(gform_ASCVD,int.lab=int.lab)

# -- for CVD
system.time(
    gform_CVD <- gfoRmula::gformula_survival(
        obs_data=datPP[datPP$ID %in% highrisk.id,], 
        id=id, time_points=time_points, time_name=time_name,
        basecovs=covnames_tf,
        covnames=covnames_tv, covtypes=covnames_tv_types,
        histories=histories_tv, histvars=histvars_tv,
        covparams=covparms_tv,
        outcome_name=outcomeY2_name, ymodel=outcomeY2_model,
        compevent_name=compeventD2_name, compevent_model=compeventD2_model,
        intvars=intvar.list, interventions=int.list, int_descript=int.desc,
        nsimul=nsimul, seed=1234, nsamples=nboot, parallel=TRUE, ncores=ncores
    )
)
print("*** Hypothetical intervenitons for CVD")
print(gform_CVD)
# coef(gform_CVD)
getRiskRatio(gform_CVD,int.lab=int.lab)

save(gform_AllCause,gform_ASCVD,gform_CVD,
     file=paste(path.res,"/pgform_Models_Results_LDLC1.8.RData",sep="")
)


load(grep("Model0",list.files(path.res),value=TRUE))
load(grep("LDLC1.8",list.files(path.res),value=TRUE))
## -- plt CIF for model 0

int.lab <- c(
  "No interention",
  "per 1 mmol/L LDL-C reduction",
  "40% reduction (Statins, HMGCR)",
  "61% reduction (PCSK9 inhibitor)")

## All-cause mortality
plt.AllCause0.CIF <- plt.CIF(
    gform.fit=gform_AllCause0,
    title="Model ASCVD risk")
plt.AllCause0.CIF
plt.AllCause.Int <- plt.Int(
    gform.fit=gform_AllCause,
    title="All-cause mortality with drug-target interventions",
    int.lab=int.lab,
    legend.idx=FALSE)
plt.AllCause.Int

## ASCVD
plt.ASCVD0.CIF <- plt.CIF(
    gform.fit=gform_ASCVD0,
    title="Model ASCVD risk")
plt.ASCVD0.CIF
plt.ASCVD.Int <- plt.Int(
    gform.fit=gform_ASCVD,
    title="ASCVD with drug-target interventions",
    int.lab=int.lab)
plt.ASCVD.Int

## CVD
plt.CVD0.CIF <- plt.CIF(
    gform.fit=gform_CVD0,
    title="Model CVD risk",
    legend.idx=FALSE)
plt.CVD0.CIF
plt.CVD.Int <- plt.Int(
    gform.fit=gform_CVD,
    title="CVD with drug-target interventionss",
    int.lab=int.lab,
    legend.idx=FALSE)
plt.CVD.Int


TTRResult <- rbind(
  cbind(getRiskRatio(gform.fit=gform_ASCVD,int.lab=int.lab),Y="ASCVD"),
  cbind(getRiskRatio(gform.fit=gform_CVD,int.lab=int.lab),Y="CVD"),
  cbind(getRiskRatio(gform.fit=gform_AllCause,int.lab=int.lab),Y="All-cause mortality")
)
TTRResult
