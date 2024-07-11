rm(list=ls())
# path.res <- this.path::this.dir()
path.res <- getwd()
suppressMessages(source(paste0(path.res,"/0_precondition_update.R")))
suppressMessages(source(paste0(path.res,"/utilities_update.R")))
suppressMessages(source(paste0(path.res,"/dat4UKBNealeLab.R")))
  
## Load variants info from UKB
selCols <- c("variant","rsid","consequence","consequence_category")
datUKBSNP <- data.table::fread(
    paste0(path.ukb,"/",grep("variant",list.files(path.ukb),value=TRUE)),
    select=selCols
)
head(datUKBSNP[1:5,])

## selection-GWAS of LDL-C from GCLC
## data set 1.
colSelected <- c("rsID","METAL_Effect","METAL_StdErr","METAL_Pvalue",
                 "POOLED_ALT_AF","ALT","REF","CHROM","POS_b37","N")
datLDLC <- data.table::fread(
    paste0(path.gclc,"/",grep("LDL_",list.files(path.gclc),value=TRUE)[1]),
    select=colSelected)
datLDLC$METAL_Pvalue <- as.numeric(datLDLC$METAL_Pvalue)
# Format the GWAS summary statistics for further analysis
datExpLDLC <- TwoSampleMR::format_data(
  as.data.frame(datLDLC),type="exposure",
  snp_col="rsID",
  beta_col="METAL_Effect",se_col="METAL_StdErr",pval_col="METAL_Pvalue",
  eaf_col="POOLED_ALT_AF",effect_allele_col="ALT",other_allele_col="REF",
  chr_col="CHROM",pos_col="POS_b37",samplesize_col="N"
)
datExpLDLC$exposure <- "GCLC trans-ancestry LDL-C"
with(datExpLDLC,table(chr.exposure))
datExpLDLC <- datExpLDLC[datExpLDLC$chr.exposure %in% datPosCtrlNew$chromosome_name,]
data.table::fwrite(
    drugTargetScreen::dat4GRAPPLE(datExpLDLC),
    file=paste(path.res,"/datExpLDLC.csv",sep="")
)

## -----------------------------------------------------------
## IVs selection strategies
## (1) genome-wide significant SNPs with LD clumped at r2<0.3
## (2) leading SNPs within each protein-coding gene locus
## (3) fine-mapping methods to select SNPs
## (4) LD clump for each gene

## (1) perform LD clump for the selection-GWAS of LDL-C
datExpLDLCClump <- drugTargetScreen::run_clump(
  dat=datExpLDLC,
  ldRef=ldRef,
  clump_kb=500,
  clump_r2=0.3,
  clump_p1=5e-08,
  threads=36
)
data.table::fwrite(
  drugTargetScreen::dat4GRAPPLE(datExpLDLCClump),
  file=paste(path.res,"/datExpLDLCClump.csv",sep="")
)

## Gene annotation for the clumped selection-GWAS of LDL-C
datExpLDLCGene <- drugTargetScreen::getAnnotatedSNPs(
  dat=datExpLDLC,
  datGRCh=datGRCh37,
  label="GCLG LDL-C"
)
data.table::setDT(datExpLDLCGene)
datExpLDLCGene <- datExpLDLCGene[,.SD[.N],by=.(SNP)]
datExpLDLCGene[
  ,`:=`(type=(ifelse(!is.na(external_gene_name) & pos>(start_position-1e5) & pos<(end_position+1e5),
                     "cis-pQTL",
                     "trans-pQTL")))
  ,by=.(SNP)]
selGenes <- unique(datExpLDLCGene$external_gene_name[datExpLDLCGene$pval<5e-06])
selGenes <- selGenes[!is.na(selGenes)]
data.table::fwrite(
  drugTargetScreen::dat4GRAPPLE(datExpLDLCGene),
  file=paste(path.res,"/datExpLDLCGene.csv",sep="")
)

## (2) obtain leading SNPs
datExpLDLCGeneLeading <- datExpLDLCGene[
  order(-pval)][,.SD[.N],by=.(external_gene_name)][
    !is.na(external_gene_name)][
      external_gene_name %in% posCtrl
    ]
with(datExpLDLCGeneLeading,addmargins(table(external_gene_name)))
names(datExpLDLCGeneLeading)[3:12] <- paste(names(datExpLDLCGeneLeading)[3:12],".exposure",sep="")
data.table::fwrite(
  drugTargetScreen::dat4GRAPPLE(datExpLDLCGeneLeading),
  file=paste(path.res,"/datExpLDLCGeneLeadning.csv",sep="")
)

## (3) Fine-mapping methods
with(datExpLDLCGene[datExpLDLCGene$external_gene_name %in% posCtrl,], table(external_gene_name))
datExpLDLCGeneSuSiE <- do.call(
  "rbind",
  lapply(posCtrl,FUN=function(i) {
    cat(paste(which(i==selGenes)," of ", length(selGenes), ": ", i, "\n", sep=""))
    datTmp <- datExpLDLC[datExpLDLC$SNP %in% datExpLDLCGene$SNP[datExpLDLCGene$external_gene_name==i],]
    n <- median(datTmp$samplesize.exposure)
    datTmpSuSiE <- drugTargetScreen::run_SuSiE(dat=datTmp,ldRef=ldRef,pval=5e-06,nSampleSize=n,threads=36,withAlleles=FALSE)
    datExpTmp <- datExpLDLC[datExpLDLC$SNP %in% datTmpSuSiE$selectedSNPs,]
    datExpTmp$external_gene_name <- i
    return(datExpTmp)
  })
)
with(datExpLDLCGeneSuSiE,addmargins(table(external_gene_name)))
## remove CETP from the positive control examples
data.table::fwrite(
  drugTargetScreen::dat4GRAPPLE(datExpLDLCGeneSuSiE[datExpLDLCGeneSuSiE$external_gene_name!="PCSK9",]),
  file=paste(path.res,"/datExpLDLCGeneSuSiE.csv",sep="")
)

## (4) LD clump for each gene
datExpLDLCGeneClump <- drugTargetScreen::run_clump(
  dat=datExpLDLC[datExpLDLC$SNP %in% datExpLDLCGene$SNP[datExpLDLCGene$external_gene_name %in% posCtrl],],
  ldRef=ldRef,
  clump_p1=1e-04,
  clump_p2=1,
  clump_kb=500,
  clump_r2=0.3,
  threads=30
)
datExpLDLCGeneClump <- merge(
  datExpLDLCGeneClump,
  datExpLDLCGene[,c("SNP","external_gene_name"),with=FALSE],
  by=c("SNP"),
  all.x=TRUE
)
with(datExpLDLCGeneClump,addmargins(table(external_gene_name)))
data.table::fwrite(
  drugTargetScreen::dat4GRAPPLE(datExpLDLCGeneClump),
  file=paste(path.res,"/datExpLDLCGeneClump.csv",sep="")
)

## ------------------------------------
## load outcome-GWAS from GARDIoGRAMC4D
datCVDOut <- as.data.frame(data.table::fread(
  paste(path.cad,"/",grep(".gz",list.files(path.cad),value=TRUE),sep="")
))
str(datCVDOut)
datCVDOut$hm_chrom <- as.integer(
  datCVDOut$hm_chrom
)
datCVDOut <- TwoSampleMR::format_data(
  datCVDOut,
  type="outcome",
  snp_col="hm_rsid",
  beta_col="hm_beta",se_col="standard_error",pval_col="p_value",
  effect_allele_col="hm_effect_allele",
  other_allele_col="hm_other_allele",
  eaf_col="hm_effect_allele_frequency",
  chr_col="hm_chrom",pos_col="hm_pos"
)
head(datCVDOut)
data.table::fwrite(
  drugTargetScreen::dat4GRAPPLE(datCVDOut),
  file=paste(path.res,"/datCVDOut.csv",sep="")
)

## harmonise exposure- and outcome-GWAS for MR analyses
## harmonise data exposure-GWAS and outcome-GWAS
# 4 versions of exposure GWAS
# (1) datMR <- datExpLDLCClump after clump with p1=5e-08, r2=0.3, window=500kb with all
# (2) datMRLeadingSNP <- datExpLDLCLeadningSNP with leading SNPs for 4 candidate genes
# (3) datMRSuSiE <- datExpLDLCSuSiE with fine-mapping via SuSiE method for each candidate gene
# (4) datMRClump <- datExpLDLCGeneClump after clump for each candidate gene

datMR <- TwoSampleMR::harmonise_data(
  datExpLDLCClump,
  datCVDOut
)
sum(datMR$mr_keep)
datMRLeadingSNP <- TwoSampleMR::harmonise_data(
  datExpLDLC[datExpLDLC$SNP %in% datExpLDLCGeneLeading$SNP,],
  datCVDOut
)
sum(datMRLeadingSNP$mr_keep)
datMRSuSiE <- TwoSampleMR::harmonise_data(
  datExpLDLCGeneSuSiE,
  datCVDOut
)
datMRSuSiE <- merge(
  datMRSuSiE,
  datExpLDLCGene[,c("SNP","type"),with=FALSE],
  by="SNP")
sum(datMRSuSiE$mr_keep)
datMRClump <- TwoSampleMR::harmonise_data(
  datExpLDLCGeneClump,
  datCVDOut
)
sum(datMRClump$mr_keep)
with(datMRClump[datMRClump$mr_keep,],
     table(external_gene_name)
)


selSNPs <- datExpLDLCGene$SNP[datExpLDLCGene$external_gene_name %in% posCtrl[posCtrl!="CETP"]]
getGeneSNPs <- function(gene) {
  datExpLDLCGene$SNP[datExpLDLCGene$external_gene_name %in% gene]
}

## Step 1. causal biomarkers screen
mainMethod <- "Inverse variance weighted"
MRResultwoCor <- drugTargetScreen::runMRIVW(datMR=datMR)[,5:9]
MRResultwoCor$exp <- "All SNPs"
MRResultwoCor[MRResultwoCor$method %in% mainMethod,]
plt.woCor <- plt.Scatter(datMR,title="IVs across the genome",legend.idx=FALSE)
plt.woCor
getMRestimates(MRResult=MRResultwoCor,mainMethod=mainMethod)


# MR results using all clumped SNPs w consideration of 4 candidate genes
MRResult4Gene <- drugTargetScreen::runMRIVW(datMR=datMR[datMR$SNP %in% selSNPs,])[,5:9]
MRResult4Gene$exp <- "APOB/LDLR/LPA"
MRResult4Gene[MRResult4Gene$method %in% mainMethod,]
plt.4Gene <- plt.Scatter(
  datMR=datMR[datMR$SNP %in% selSNPs,],
  title="IVs within APOB/LDLR/LPA loci")
plt.4Gene
getMRestimates(MRResult=MRResult4Gene,mainMethod=mainMethod)

save(datExpLDLCGene,datMR,datMRLeadingSNP,datMRSuSiE,datMRClump,
     MRResultwoCor,plt.woCor,
     MRResult4Gene,plt.4Gene,
     file=paste0(path.C,"/RealExample_1_causal_biomarkers_screen_20240709.RData"))
