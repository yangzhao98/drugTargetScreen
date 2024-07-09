#' @title rename TwoSampleMR::format_data() with snps selection at \code{pval}
#'
#' @param dat data frame obtained from TwoSampleMR::format_data()
#' @param label label the data frame
#' @param pval p-value threshod for selecting significant snps
#'
#' @export
rmName4formatDat <- function(dat,label,pval) {
  SNP <- chr <- pos <- NULL
  names(dat) <- gsub(".exposure","",names(dat))
  names(dat) <- gsub(".outcome","",names(dat))
  dat <- subset(dat,select=c(SNP,chr,pos,pval))
  dat$chr <- with(dat,as.numeric(as.character(chr)))
  dat <- dat[!is.na(dat$chr) & dat$pval < pval,]
  names(dat) <- c(names(dat)[1:3],label)
  return(dat)
}


#' @title merge the data set based on SNP
#'
#' @param datX data frame X. Usually, from \code{rmName4formatDat()}
#' @param datY data frame Y. Usually, from \code{rmName4FormatDat()}
#'
#' @export
datMerge <- function(datX,datY) {
  SNP <- chr <- pos <- NULL
  snpList <- dplyr::intersect(datX$SNP,datY$SNP)
  snpListX<- dplyr::setdiff(datX$SNP,datY$SNP)
  snpListY<- dplyr::setdiff(datY$SNP,datX$SNP)

  dat <- merge(datX[datX$SNP %in% snpList,],
               datY[datY$SNP %in% snpList, ],
               by = c("SNP","chr","pos"))
  dat <- plyr::rbind.fill(
    dat,
    datX[datX$SNP %in% snpListX, ],
    datY[datY$SNP %in% snpListY, ])
  return(dat)
}


#' @title annotate the reName4formatDat with gene info
#'
#' @param rmNameformatDat processed TwoSampleMR::format_data()
#' @param datGRCh reference data set with annotation info
#'
#' @export
annotateGeneName <- function(rmNameformatDat,datGRCh) {
  chr <- pos <- chromosome_name <- start_position <- end_position <- NULL
  dplyr::left_join(
    rmNameformatDat,datGRCh,
    dplyr::join_by(chr==chromosome_name,
                   pos>=start_position,
                   pos<=end_position))
}



#' @title Get annotated SNPs info
#'
#' @param dat data frame obtained from \code{TwoSampleMR::format_data()} for exposure GWAS
#' @param datGRCh reference panel for gene annotation in either GRCh37 or GRCh38
#' @param label label the trait or phenotype info
#'
#' @export
getAnnotatedSNPs <- function(dat,datGRCh,label) {
  ## formated exposure GWAS data using TwoSampleMR::format_data()
  dat <- dat[,!names(dat) %in% c("pval_origin.exposure","id.exposure","exposure")]
  dat$trait <- label
  ## annotation
  datGene <- drugTargetScreen::annotateGeneName(
    drugTargetScreen::rmName4formatDat(dat=dat,label=label,pval=1),
    datGRCh=datGRCh
  )
  datGene <- datGene[,c(1,5:10)]
  ## merge data
  datGene <- merge(dat, datGene,by="SNP")
  names(datGene) <- gsub(".exposure","",names(datGene))
  # names(datGene)[2:3] <- c("chr_GRCh37","pos_GRCh37")
  datGene <- datGene[order(datGene$chr,datGene$pos),]
  return(datGene)
}




#' @title Implement the fine-mapping analysis using SuSiE method
#'
#' @param dat data frame from \code{TwoSampleMR::format_data()} for exposure GWAS
#' @param ldRef location of 1KG reference panel
#' @param pval pvalue threshold for selecting SNPs
#' @param nSampleSize sample size of the exposure GWAS
#' @param threads number of threads used for plink
#' @param withAlleles indicator of whether LD matrix includes effect and other alleles
#'
#' @export
#'
run_SuSiE <- function(dat,ldRef,pval,nSampleSize,threads,withAlleles=TRUE) {
  # get LD matrix for PSCK9 gene
  datTmp <- getLDmatrix(dat=dat,ldRef=ldRef,pval=pval,
                        threads=threads,withAlleles=withAlleles)
  # calc the z-score
  datTmp$dat$zscore <- with(datTmp$dat,beta.exposure/se.exposure)
  # run SuSiE
  fittedSuSiE <- tryCatch(
    with(
      datTmp$dat,
      susieR::susie_rss(z=zscore,
                        R=datTmp$LDmatrix,
                        L=10,
                        n=nSampleSize,
                        estimate_residual_variance=TRUE)),
    error=function(e) {
      with(
        datTmp$dat,
        susieR::susie_rss(z=zscore,
                          R=datTmp$LDmatrix,
                          L=10,
                          n=nSampleSize,
                          estimate_residual_variance=FALSE))
    }
  )

  calPIP <- summary(fittedSuSiE)$vars
  selectedSNPidx <- sort(calPIP$variable[calPIP$cs>0])
  # plot the regional plot
  pltSuSiE <- susieR::susie_plot(
    fittedSuSiE,y="PIP",
    b=datTmp$dat$beta.exposure,
    xlab="PCSK9")
  # plot the LD plot
  pltLD <- gaston::LD.plot(
    datTmp$LDmatrix,
    snp.positions=datTmp$dat$pos.exposure,
    polygon.par = list(border = NA),
    write.ld = NULL)

  return(list(dat=datTmp$dat,LDmatrix=datTmp$LDmatrix,
              selectedSNPs=datTmp$dat$SNP[selectedSNPidx],
              pltSuSiE=pltSuSiE,
              pltLD=pltLD))

}



#' @title remove ".exposure" and ".outcome" from a data frame
#' @param dat GWAS summary statistics from \code{TwoSampleMR::format_data()}
#' @keywords internal
#' @export
rnDat <- function(dat) {
  names(dat) <- gsub(".exposure","",names(dat))
  names(dat) <- gsub(".outcome","",names(dat))
  return(dat)
}


#' @title calculate number of the overlapped SNPs between two traits given a pval threshold
#' @param dat1 GWAS summary statistics from \code{TwoSampleMR::format_data()} for trait 1
#' @param dat2 GWAS summary statistics from \code{TwoSampleMR::format_data()} for trait 2
#' @param pval pvalue threshod for selecting the overlapped SNPs
#' @keywords internal
#' @export
checkSNPsOverlap <- function(dat1,dat2,pval) {
  dat1 <- rnDat(dat1)
  dat2 <- rnDat(dat2)
  snps <- intersect(dat1$SNP[dat1$pval<pval],dat2$SNP[dat2$pval<pval])
  return(list(nSNPs=length(snps),SNPs=snps))
}

#' @title get the overlapped SNPs between two traits for colocalization analyses
#' @param dat1 GWAS summary statistics from \code{TwoSampleMR::format_data()} for trait1
#' @param trait1 trait name for dat1
#' @param dat2 GWAS summary statistics from \code{TwoSAmpleMR::format_data()} for trait2
#' @param trait2 trait name for dat2
#' @param ldRef 1KG reference panel in plink format
#' @param pval pvalue threshold for selecting SNPs with \code{pvalue<pval}
#' @keywords internal
#' @export
getOverlapSNPs <- function(dat1,trait1,dat2,trait2,ldRef,pval) {
  ## Strategy for getting overlapped SNPs
  # (1) preliminary check the overlapped SNPs between trait1 and trait2
  # (2) check the overlapped SNPs in 1KG reference and get LDmatrix
  # (3) select the overlapped SNPs in trait1, trait2, and 1KG reference panel
  dat1 <- rnDat(dat1)
  dat2 <- rnDat(dat2)
  SNPs <- checkSNPsOverlap(dat1=dat1,dat2=dat2,pval=pval)
  ## Check SNPs in 1KG reference
  dat1 <- drugTargetScreen::getLDmatrix(dat1[dat1$SNP %in% SNPs$SNPs,],ldRef=ldRef,pval=pval)
  dat2 <- drugTargetScreen::getLDmatrix(dat2[dat2$SNP %in% SNPs$SNPs,],ldRef=ldRef,pval=pval)
  snps <- intersect(dat1$dat$SNP,dat2$dat$SNP)
  ## Get the overlapp SNPs
  dat1$dat <- dat1$dat[dat1$dat$SNP %in% snps,]
  rownames(dat1$LDmatrix) <- unlist(lapply(rownames(dat1$LDmatrix),FUN=function(x) strsplit(x,"_")[[1]][1]))
  colnames(dat1$LDmatrix) <- unlist(lapply(colnames(dat1$LDmatrix),FUN=function(x) strsplit(x,"_")[[1]][1]))
  dat1$LDmatrix <- dat1$LDmatrix[rownames(dat1$LDmatrix) %in% snps,colnames(dat1$LDmatrix) %in% snps]
  dat2$dat <- dat2$dat[dat2$dat$SNP %in% snps,]
  rownames(dat2$LDmatrix) <- unlist(lapply(rownames(dat2$LDmatrix),FUN=function(x) strsplit(x,"_")[[1]][1]))
  colnames(dat2$LDmatrix) <- unlist(lapply(colnames(dat2$LDmatrix),FUN=function(x) strsplit(x,"_")[[1]][1]))
  dat2$LDmatrix <- dat2$LDmatrix[rownames(dat2$LDmatrix) %in% snps,colnames(dat2$LDmatrix) %in% snps]
  ## Get the overlapp LDmatrix
  return(stats::setNames(list(dat1,dat2),c(trait1,trait2)))
}


#' @title prepare data for implementing MR-GRAPPLE
#' @param dat formated exposue- or outcome-GWAS using TwoSampleMR::format_data()
#' @export
dat4GRAPPLE <- function(dat) {
  names(dat) <- gsub(".exposure","",names(dat))
  names(dat) <- gsub(".outcome","",names(dat))
  chkCols <- c("SNP","effect_allele","other_allele","beta","se","pval")
  if (length(match.arg(names(dat),chkCols,several.ok=TRUE))!=6) {
    stop(paste("Please check the input data with at least \n",chkCols,sep=""))
  }
  return(dat)
}

#' @title data processing for UKB-Neale Labs GWAS summary statistics
#' @param UKBGWAS.file path of the GWAS summary statistics
#' @param datUKBSNP processed variants info obtained from variants.tsv.gz with \code{c("variant","rsid","consequence","consequence_category")}
#' @param type the type of exposure- or outcome- GWAS processed using \code{TwoSampleMR::format_data()}. By default, \code{type="outcome"}
#' @param binaryTrait an indicator of a binary trait or not. By default, \code{binaryTrait=TRUE}.
#' @import data.table
#' @export
dat4UKBNealeLab <- function(UKBGWAS.file,datUKBSNP,type="outcome",binaryTrait=TRUE) {
  ## Variants info
  # selCols <- c("variant","rsid","consequence","consequence_category")
  # datUKBSNP <- data.table::fread(UKBSNP.file,select=selCols)
  ## GWAS info
  datUKBGWAS <- data.table::fread(UKBGWAS.file)
  datUKBGWAS <- datUKBGWAS[datUKBGWAS$low_confidence_variant==FALSE,]
  ## log(OR) transformation from BOLT-LMM binary trait analysis
  ## ref: https://data.bris.ac.uk/data/dataset/pnoat8cxo0u52p6ynfaekeigi
  logORTransformation <- function(beta,nCase,nTot) {
    u <- nCase/nTot
    beta/(u*(1-u))
  }
  ##
  data.table::setDT(datUKBGWAS);data.table::setkey(datUKBGWAS,"variant")
  data.table::setDT(datUKBSNP); data.table::setkey(datUKBSNP,"variant")
  datUKBGWAS <- datUKBGWAS[datUKBSNP,nomatch=0]
  #system.time(datUKBGWAS <- merge(datUKBGWAS,datUKBSNP,by="variant",all.x=TRUE))
  a <- strsplit(datUKBGWAS$variant,":")
  datUKBGWAS$chr <- sapply(a,"[[",1)
  datUKBGWAS$pos <- sapply(a,"[[",2)
  datUKBGWAS$ea <- sapply(a,"[[",3)
  datUKBGWAS$ra <- sapply(a,"[[",4)
  if (binaryTrait) {
    datUKBGWAS[
      ,`:=`(nCase=ceiling(expected_case_minor_AC/(2*minor_AF)))
      ,by=.(variant)]
    datUKBGWAS[
      ,`:=`(effect_allele=ifelse(ea==minor_allele,ea,ra),
            eaf=ifelse(ea==minor_allele,minor_AF,1-minor_AF),
            other_allele=ifelse(ea==minor_allele,ra,ea),
            betaNew=logORTransformation(beta,nCase,n_complete_samples),
            seNew=logORTransformation(se,nCase,n_complete_samples))
      ,by=.(variant)][,unique("variant")]
    datUKBGWAS <- TwoSampleMR::format_data(
      dat=as.data.frame(datUKBGWAS),
      type=type,
      snp_col="rsid",beta_col="betaNew",se_col="seNew",pval_col="pval",
      effect_allele_col="effect_allele",
      other_allele_col="other_allele",
      eaf_col="eaf",
      chr_col="chr",pos_col="pos",samplesize_col="n_complete_samples")
  } else {
    datUKBGWAS[
      ,`:=`(effect_allele=ifelse(ea==minor_allele,ea,ra),
            eaf=ifelse(ea==minor_allele,minor_AF,1-minor_AF),
            other_allele=ifelse(ea==minor_allele,ra,ea))
      ,by=.(variant)][,unique("variant")]
    datUKBGWAS <- TwoSampleMR::format_data(
      dat=as.data.frame(datUKBGWAS),
      type=type,
      snp_col="rsid",beta_col="beta",se_col="se",pval_col="pval",
      effect_allele_col="effect_allele",
      other_allele_col="other_allele",
      eaf_col="eaf",
      chr_col="chr",pos_col="pos",samplesize_col="n_complete_samples")
  }
  datUKBGWAS <- datUKBGWAS[!is.na(datUKBGWAS$SNP),]
  return(datUKBGWAS)
}
