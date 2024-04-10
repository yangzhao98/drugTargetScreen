#' @title plot the regional association plot for exp- and out-GWAS
#' @param geneInfo dataset contains the gene info that is generally obtained from \code{datGRCh37} or \code{datGRCh38}.
#' @param datOut out-GWAS summary statistics obtained from \code{TwoSampleMR::format_data()}
#' @param outName name of the outcome trait
#' @param datExp exp-GWAS summary statistics obtained from \code{TwoSampleMR::format_data()}
#' @param expName name of the exposure trait
#' @param ldRef 1KG reference panel in plink format
#' @export
plt.region <- function(geneInfo,datOut,outName,datExp,expName,ldRef) {
  ## ----------------------------------------------------------------------
  ## Strategies for plot region
  ## (1) select the overlapped SNPs
  ## (2) calculate the correlation between SNPs using ldRef
  ## (3) prepare the data with head c(marker,chr,pos,pvalue_1,pvalue_2,...)
  ## (4) provide the trait names
  ## ----------------------------
  ## gene info
  chrpos <- with(geneInfo,c(chromosome_name,start_position,end_position))
  ## get overlapped dataset
  datTmp <- drugTargetScreen::getOverlapSNPs(
    dat1=datExp[datExp$chr.exposure==chrpos[1] &
                  datExp$pos.exposure>=chrpos[2]-5e4 &
                  datExp$pos.exposure<=chrpos[3]+5e4,],
    trait1=expName,
    dat2=datOut,trait2=outName,
    ldRef=ldRef,pval=1)
  ## prepare data for stacked coloc plots
  ## Note: the same data format for CMplot
  datColoc <- merge(datTmp[[expName]]$dat[,c("SNP","chr","pos","pval")],
                    datTmp[[outName]]$dat[,c("SNP","pval")],
                    by="SNP")
  names(datColoc) <- c("marker","chr","pos",paste("pvalue_",1:(ncol(datColoc)-3),sep=""))
  datColoc$chr <- as.integer(datColoc$chr)
  datColoc$pos <- as.integer(datColoc$pos)
  datColoc <- datColoc[order(datColoc$chr,datColoc$pos),]
  return(geni.plots::fig_region_stack(
    data=datColoc,
    traits=c(expName,outName),
    corr=datTmp[[expName]]$LDmatrix,
    build=37,
    title_center=TRUE
  ))
}


#' @title plot the regional association plot for exp-, eQTL- and out-GWAS
#' @param geneInfo dataset contains the gene info that is generally obtained from \code{datGRCh37} or \code{datGRCh38}.
#' @param datOut out-GWAS summary statistics obtained from \code{TwoSampleMR::format_data()}
#' @param outName name of the outcome trait
#' @param datQTL QTL-GWAS summary statistics obtained from \code{TwoSampleMR::format_data()} or with the same data format
#' @param qtlName name of the eQTL trait
#' @param datExp exp-GWAS summary statistics obtained from \code{TwoSampleMR::format_data()}
#' @param expName name of the exposure trait
#' @param ldRef 1KG reference panel in plink format
#' @export
plt.stack.region <- function(geneInfo,datOut,outName,
                             datQTL,qtlName,datExp,expName,
                             ldRef) {
  ## gene info
  chrpos <- with(geneInfo,c(chromosome_name,start_position,end_position))
  ## get overlapped SNPs
  # (1) get overlapped SNPs between exp-GWAS and out-GWAS
  datExpOut <- drugTargetScreen::getOverlapSNPs(
    dat1=datExp[datExp$chr.exposure==chrpos[1] &
                  datExp$pos.exposure>=chrpos[2]-5e4 &
                  datExp$pos.exposure<=chrpos[3]+5e4,],
    trait1=expName,
    dat2=datOut,trait2=outName,
    ldRef=ldRef,pval=1)
  # (2) get overlapped SNPs between exp-GWAS and eQTL-GWAS
  datExpOuteQTL <- drugTargetScreen::getOverlapSNPs(
    dat1=datExpOut[[expName]]$dat,trait1=expName,
    dat2=datQTL,trait2=qtlName,
    ldRef=ldRef,pval=1
  )
  snp <- datExpOuteQTL[[expName]]$dat$SNP
  ldMatrix <- datExpOuteQTL[[expName]]$LDmatrix
  datColoc <- merge(merge(
    datExpOuteQTL[[expName]]$dat[,c("SNP","chr","pos","pval")],
    datExpOuteQTL[[qtlName]]$dat[,c("SNP","pval")],by="SNP"),
    datExpOut[[outName]]$dat[,c("SNP","pval")],by="SNP"
  )
  names(datColoc) <- c("marker","chr","pos",paste("pvalue_",1:(ncol(datColoc)-3),sep=""))
  # head(datColoc)
  snpList <- row.names(ldMatrix)
  idx <- unlist(lapply(snpList, function(i) which(i==datColoc$marker)))
  datColoc <- datColoc[idx,]
  return(geni.plots::fig_region_stack(
    data=datColoc,
    traits=c(expName,qtlName,outName),
    corr=ldMatrix,
    build=37,
    title_center=TRUE
  ))
}

