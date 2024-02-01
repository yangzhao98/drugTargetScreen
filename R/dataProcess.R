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
