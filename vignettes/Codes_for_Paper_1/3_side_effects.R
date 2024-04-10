# devtools::install_github("yangzhao98/drugTargetScreen",force=TRUE)

datUKBTab <- readxl::read_xlsx(
    paste(path.PheMR,"/",
          grep(".xlsx",list.files(path.PheMR),value=TRUE),
          sep=""
    ),sheet=2
)
data.table::setDT(datUKBTab)
datUKBTab <- datUKBTab[grepl("Diagnoses - main ICD10",datUKBTab$`Phenotype Description`),]
datUKBTab <- datUKBTab[grepl("both_sexes",Sex)]
datUKBTab <- datUKBTab[,`:=`(sidePhe=trimws(gsub(paste("Diagnoses - main ICD10: ",`Phenotype Code`,sep=""),"",
                                          `Phenotype Description`)),which="both")
                       ,by=.(`Phenotype Code`)
    ][,c("Phenotype Code","sidePhe","Phenotype Description"),with=FALSE
]
datPheMR.files <- list.files(paste(path.PheMR,"/UKB_ICD10_MainDiagnosis",sep=""))
UKBSNP.file <- paste(path.PheMR,"/",grep("variant",list.files(path.PheMR),value=TRUE),sep="")
UKBGWAS.file <- paste(path.PheMR,"/UKB_ICD10_MainDiagnosis/",datPheMR.files[1],sep="")

# datUKBOut <- dat4UKBNealeLab(UKBGWAS.file=UKBGWAS.file,datUKBSNP=datUKBSNP)
MRResultPheMR <- do.call(
    "rbind",
    lapply(datPheMR.files,FUN=function(i) {
        cat(paste(which(datPheMR.files==i), " of ", length(datPheMR.files), ": ",i,"\n",
                  as.character(Sys.time()),"\n",sep=""))
        UKBGWAS.file <- paste(path.PheMR,"/UKB_ICD10_MainDiagnosis/",i,sep="")
        datUKBOut <- dat4UKBNealeLab(UKBGWAS.file=UKBGWAS.file,datUKBSNP=datUKBSNP)
        MRRes <- NULL
        for (j in posCtrl) {
            datExpTmp <- datExpLDLCGeneSuSiE[datExpLDLCGeneSuSiE$external_gene_name==j,]
            datMRTmp <- TwoSampleMR::harmonise_data(datExpTmp,datUKBOut)
            if(nrow(datMRTmp)>0 & sum(datMRTmp$mr_keep)>0) {
                MRResultTmp <- drugTargetScreen::runMRIVW(dat=datMRTmp)
                MRResultTmp$external_gene_name =j
                MRRes <- rbind(MRRes,MRResultTmp)
            } else {
                MRResultTmp <- NULL
                MRRes <- rbind(MRRes,MRResultTmp)
            }
        }
        if (nrow(MRRes)>0) MRRes$UKBOut <- i else MRRes <- NULL
        return(MRRes)
    })
)

## clear the results
load(paste(path.res,"/MRResultPheMR.RData",sep=""))
data.table::setDT(MRResultPheMR)
MRResultPheMR <- MRResultPheMR[
  ,`:=`(ukbID=strsplit(UKBOut,"\\.",)[[1]][1])
  ,by=.(UKBOut,method)]
MRResultPheMR <- merge(datUKBTab,MRResultPheMR,by.y="ukbID",by.x="Phenotype Code")
FDR.level <- 0.05
MRResultPheMR[,`:=`(log10p=log10(as.numeric(pval)),OR=exp(b))]
MRResultPheMR[,`:=`(qval=qvalue::qvalue(p=as.numeric(pval),fdr.level=FDR.level)$qvalues)
              ,by=.(external_gene_name,method)]
MRResultPheMR[,`:=`(idx=qval<FDR.level,
                    dir=ifelse(b<0,"Protective","Hazard"))]
## plt volcano plot
plt.4Gene.Volcano <- lapply(posCtrl,function(i) plt.Volcano(mainMethod,gene=i))


