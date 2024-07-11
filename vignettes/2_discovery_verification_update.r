## conduct verious MR analyses
# (1) MR with IVs selected from LD clumped SNPs for each gene
# (2) MR with IVs selected from leading SNP for each gene
# (3) MR with IVs selected from fine-mapping method of SuSiE
# (4) MR with cluster method
# (5) MR with GRAPPLE in a three-sample setting
# (6) MR with IVs identified cis-pQTLs


## ---
# (1) MR results using clumpled SNPs obtained from each candidate gene
MRResultClump <- do.call(
  "rbind",
  lapply(posCtrl,FUN=function(i){
    datMRTmp <- datMRClump[datMRClump$external_gene_name %in% i,]
    if (nrow(datMRTmp)>0) {
      datMRRes <- drugTargetScreen::runMRIVW(datMR=datMRTmp)[,5:9]
      datMRRes$exp <- paste("Clump SNPs of ", i, sep="")
    } else {
      datMRRes <- NULL
    }
    datMRRes
  })
)
MRResultClump[MRResultClump$method %in% mainMethod,]
plt.Clump <- lapply(posCtrl,FUN=function(x) {
  plt.Scatter(datMRClump[datMRClump$external_gene_name %in% x,],
              title=paste("IVs within ",x, " locus"))
}
)
plt.Clump[[1]]
plt.Clump[[2]]
plt.Clump[[3]]
plt.Clump[[4]]


## ---
# (2) MR results using leading SNP for each gene
MRResultLeading <- do.call(
  "rbind",
  lapply(posCtrl,FUN=function(i) {
    datMRTmp <- datMRLeadingSNP[
      datMRLeadingSNP$SNP %in% datExpLDLCGeneLeading$SNP[datExpLDLCGeneLeading$external_gene_name==i],
    ]
    if (nrow(datMRTmp)>0) {
      datMRRes <- drugTargetScreen::runMRIVW(datMR=datMRTmp)[,5:9]
      datMRRes$exp <- paste("Leading SNP of ", i, sep="")
    } else {
      datMRRes <- NULL
    }
    datMRRes
  })
)
MRResultLeading
plt.Leading <- plt.Scatter(datMRLeadingSNP,title="Leading SNPs")
plt.Leading

## ---
# (3) MR results suing SNPs obtained from fine-mapping using SuSiE
MRResultSuSiE <- do.call(
  "rbind",
  lapply(posCtrl,FUN=function(i){
    datMRTmp <- datMRSuSiE[datMRSuSiE$external_gene_name %in% i,]
    if (nrow(datMRTmp)>0) {
      datMRRes <- drugTargetScreen::runMRIVW(datMR=datMRTmp)[,5:9]
      datMRRes$exp <- paste("Fine-mapping SNPs of ", i, sep="")
    } else {
      datMRRes <- NULL
    }
    datMRRes
  })
)
getMRestimates(MRResult=MRResultSuSiE,mainMethod=mainMethod)
MRResultSuSiE[MRResultSuSiE$method %in% mainMethod,]
plt.SuSiE <- lapply(posCtrl,FUN=function(x) {
  plt.Scatter(datMRSuSiE[datMRSuSiE$external_gene_name %in% x,],
              title=paste("IVs (SuSiE) within ",x, " locus"))
}
)
plt.SuSiE[[1]]
plt.SuSiE[[2]]
plt.SuSiE[[3]]
plt.SuSiE[[4]]


### ---
# (4) MRClust method
# MRResultAllSNPs <- drugTargetScreen::runMRClust(datMR=datMR)
# plt.Clust <- MRResultAllSNPs$pltScatter
# plt.Clust
# MRResultAllSNPs$pltScatterProb08.2IVs
# MRResultSuSiEClust <- runMRClust(datMR=datMRSuSiE)
# MRResultClumpClust <- runMRClust(datMR=datMRClump)

## --------------------------------------------------------------------------------------------------------
## -----------------------------Step 2. Run GRAPPLE with three samples MR----------------------------------
selGWAS <- paste(path.res,"/datExpLDLCGeneSuSiE.csv",sep="")
expGWAS <- paste(path.res,"/datExpUKB.csv",sep="")
outGWAS <- paste(path.res,"/datCVDOut.csv",sep="")
# datExpUKB <- data.table::fread(expGWAS)
MRGRAPPLE <- drugTargetScreen::runGRAPPLE(
  selGWAS.file=selGWAS,
  expGWAS.file=expGWAS,
  outGWAS.file=outGWAS,
  ldRef=ldRef,
  pval=5e-06,
  clump_r2=0.8,
  outPath=path.res
)
with(MRGRAPPLE$datGRAPPLE.ProteinCodingGene,
     table(GENCODE_name,Mode1_marker)
)
mode1 <- MRGRAPPLE$datGRAPPLE.ProteinCodingGene$rsID[MRGRAPPLE$datGRAPPLE.ProteinCodingGene$Mode1_marker]
mode1 <- mode1[!is.na(mode1)]
mode2 <- MRGRAPPLE$datGRAPPLE.ProteinCodingGene$rsID[MRGRAPPLE$datGRAPPLE.ProteinCodingGene$Mode2_marker]
mode2 <- mode2[!is.na(mode2)]
datMRGRAPPLE <- data.frame(
  SNP=c(mode1,mode2),
  Model=c(rep(1,length(mode1)),rep(2,length(mode2)))
)
datMRGRAPPLE <- merge(
  datMRGRAPPLE,
  MRGRAPPLE$datGRAPPLE.ProteinCodingGene[,c("rsID","GENCODE_name")],
  by.x="SNP",by.y="rsID",all.x=TRUE
)
datMRGRAPPLE$idx <- with(datMRGRAPPLE,paste(GENCODE_name,"-",Model,sep=""))
idx <- unique(datMRGRAPPLE$idx)
datMRSuSiE[datMRSuSiE$SNP %in% datMRGRAPPLE$SNP,]
MRGRAPPLEMode <- do.call(
  "rbind",
  lapply(posCtrl[!posCtrl %in% c("PCSK9")],FUN=function(i) {
    datMRTmp <- datMRSuSiE[datMRSuSiE$SNP %in% datMRGRAPPLE$SNP[datMRGRAPPLE$GENCODE_name==i],]
    if (nrow(datMRTmp)>0) {
      datMRRes <- drugTargetScreen::runMRIVW(datMR=datMRTmp)[,5:9]
      datMRRes$exp <- paste("MR-GRAPPLE: ", i, sep="")
    } else {
      datMRRes <- NULL
    }
    datMRRes
  })
)
# mainMethod <- c("Inverse variance weighted")
# MRGRAPPLEMode[MRGRAPPLEMode$method %in% mainMethod,]
plt.GRAPPLE <- lapply(1:2,FUN=function(x) {
  plt.Scatter(datMRSuSiE[datMRSuSiE$SNP %in% datMRGRAPPLE$SNP[datMRGRAPPLE$Model==x],],
              title=paste("IVs included in Model ",x,sep=""))
}
)
plt.GRAPPLE[[1]]
plt.GRAPPLE[[2]]
plt.GRAPPLE.profile <- MRGRAPPLE$pltGRAPPLE +
  ggplot2::theme(axis.ticks.y=element_line(colour="black"),
                 axis.text.y=element_text(colour="black"),
                 plot.title=element_text(size=14,hjust=0.5),
                 axis.text=element_text(size=14),
                 axis.title=element_text(size=14)) +
  ggplot2::scale_x_continuous(name=expression(~beta))
plt.GRAPPLE.profile

## cis- and trans-MR analyses
dim(datMRClump)
datMRClump <- merge(
  datMRClump,
  datExpLDLCGene[,c("SNP","type"),with=FALSE],
  by="SNP")

MRResultCis <- do.call(
  "rbind",
  lapply(posCtrl,FUN=function(i){
    datMRTmp <- datMRClump[datMRClump$external_gene_name %in% i,]
    if (nrow(datMRTmp)>0) {
      datMRRes <- drugTargetScreen::runMRIVW(datMR=datMRTmp)[,5:9]
      datMRRes$exp <- paste("Fine-mapping SNPs of ", i, sep="")
    } else {
      datMRRes <- NULL
    }
    datMRRes
  })
)
MRResultCis[MRResultCis$method %in% mainMethod,]
plt.Cis <- lapply(posCtrl,FUN=function(x) {
  plt.Scatter(datMRClump[datMRClump$external_gene_name %in% x,],
              title=paste("Cis-pQTLs within ",x, " locus"))
}
)
plt.Cis[[1]]
plt.Cis[[2]]
plt.Cis[[3]]
plt.Cis[[4]]


## Plot forest
datMRForest <- function(MRResult,label) {
  selColumns <- c("method","nsnp","b","se","pval","exp")
  data.table::setDT(MRResult)
  MRResult <- MRResult[,selColumns,with=FALSE]
  MRResult[,`:=`(label=label)]
  return(as.data.frame(MRResult))
}
dat4Forest <- rbind(
  datMRForest(MRResult4Gene,label="All IVs"),
  datMRForest(MRResultSuSiE,label="SuSiE"),
  datMRForest(MRGRAPPLEMode,label="GRAPPLE"),
  datMRForest(MRResultCis,label="Cis-pQTLs")
)
data.table::setDT(dat4Forest)
dat4Forest <- dat4Forest[!label %in% c("All IVs","CETP","PCSK9")]
# dat4Forest <- dat4Forest[!Gene %in% c("All IVs","CETP")]
dat4Forest[,`:=`(Gene=ifelse(label!="All IVs" & grepl("APOB",exp),"APOB",
                             ifelse(label!="All IVs" & grepl("LDLR",exp),"LDLR",
                                    ifelse(label!="All IVs" & grepl("LPA",exp),"LPA",
                                           ifelse(label!="All IVs" & grepl("PCSK9",exp),"PCSK9","All")))))
           ,by=.(exp)]
dat4Forest[,`:=`(or=exp(b),lci=exp(b-1.96*se),uci=exp(b+1.96*se))]
mainMethod <- c("Inverse variance weighted","Wald ratio")
dat4ForestIVW <- dat4Forest[method %in% mainMethod]
library(ggplot2)
plt.Forest <- ggplot2::ggplot(
  dat4ForestIVW[!Gene %in% c("PCSK9")],
  ggplot2::aes(y=or,ymin=lci,ymax=uci,x=Gene,col=label)) +
  ggplot2::geom_pointrange(position=position_dodge(width=0.6)) +
  ggplot2::theme_classic() +
  ggplot2::geom_hline(aes(yintercept=1),linetype=2) +
  ggplot2::scale_y_continuous(transform = "log",breaks=exp(log(c(0,1,2.5,5,7.5,10,15)))) +
  ggplot2::coord_flip() +
  ggplot2::labs(y="Odds ratio & 95% confidence intervals",
                x="Drug-targted genes") +
  ggplot2::theme(legend.position=c(0.8,0.3),
                 legend.justification='top',
                 legend.box.just="right",
                 legend.direction="vertical",
                 legend.title=element_blank(),
                 axis.line = element_line(linetype = "solid"),
                 panel.background=element_blank(),
                 panel.border=element_blank(),
                 panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 plot.background=element_blank(),
                 plot.title=element_text(size=14,hjust=0.5),
                 axis.text=element_text(size=14),
                 axis.title=element_text(size=14)) +
  ggsci::scale_color_jama() +
  ggtitle("Mimicked effects of per SD increase in LDL-C on CAD")
plt.Forest

## ----------------------------------------------------------------------------------------
## -- Additional layers
## -- Colocalization with eQTLs
datPosCtrlNew
selCols.eQTL <- c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","NrSamples","GeneSymbol")
dateQTLCis <- data.table::fread(
  paste(path.gtexv8,"/",grep("2019-12-11",list.files(path.gtexv8),value=TRUE),sep=""),
  select=selCols.eQTL
)
posCtrl <- datPosCtrlNew$external_gene_name
# head(dateQTLNew)
dateQTL <- dateQTLCis
names(dateQTL) <- c(
  "pval.exposure","SNP","chr.exposure","pos.exposure",
  "effect_allele.exposure","other_allele.exposure","samplesize.exposure","external_gene_name"
)
dateQTL <- dateQTL[!duplicated(dateQTL$SNP),]

datExpLDLC$chr.exposure <- as.integer(datExpLDLC$chr.exposure)
datExpLDLC$pos.exposure <- as.integer(datExpLDLC$pos.exposure)

## original GCLC
posCtrlNew <- c("APOB","LDLR","LPA")
plt.4Gene.CAD.runColoc <- lapply(
  posCtrlNew,FUN=function(i) {
    cat(paste0("Gene: ",i,"\n"))
    runColoc(
      geneInfo=datPosCtrlNew[external_gene_name==i],
      datOut=datCVDOut,outName="Coronary artery disease",
      datQTL=dateQTL,qtlName="eQTLGene phase 1",
      datExp=datExpLDLC,expName=paste("Mimicked effects of per SD increase in LDL-C by ",i,sep=""),
      distanceGene=1e4)
  }
)
plt.4Gene.CAD.runColoc[[1]]
plt.4Gene.CAD.runColoc[[2]]
plt.4Gene.CAD.runColoc[[3]]


save(MRResultClump,plt.Clump,
     MRResultLeading,plt.Leading,
     MRResultSuSiE,plt.SuSiE,
     MRGRAPPLE,MRGRAPPLEMode,plt.GRAPPLE.profile,
     MRResultCis,plt.Cis,
     datMRForest,dat4Forest,plt.Forest,
     plt.4Gene.CAD.runColoc,
     file=paste0(path.C,"/RealExample_2_discovery_verfication_20240709.RData"))
