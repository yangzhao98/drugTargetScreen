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


#' @title run ld clump locally
#'
#' @param dat data frame from \code{TwoSampleMR::format_data()} for exposure GWAS
#' @param ldRef 1KG reference panel in plink format
#' @param clump_p1 p-value threshold for the index SNP. By default, \code{clump_p1=1}
#' @param clump_p2 p-value threshold for the secondary SNP. By default, \code{clump_p2=1}
#' @param clump_kb clumping window. By default, \code{clump_kb=10000}
#' @param clump_r2 clumping r2 cutoff. By default, \code{clump_r2=0.001}
#' @param threads number of threads used for plink
#'
#' @export
run_clump <- function(dat,
                      ldRef,
                      clump_kb = 10000,
                      clump_r2 = 0.001,
                      clump_p1 = 5e-08,
                      clump_p2 = 1,
                      threads = 36) {
  SNP <- pval.exposure <- NULL
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")
  fn <- tempfile()
  utils::write.table(data.frame(SNP=dat[["SNP"]],P=dat[["pval.exposure"]]),
                     file = fn, row.names = F, col.names = T, quote = F)
  fun2 <- paste0(shQuote(plinkbinr::get_plink_exe(), type = shell),
                 " --bfile ", shQuote(ldRef, type = shell),
                 " --clump ", shQuote(fn, type = shell),
                 " --clump-p1 ", clump_p1,
                 " --clump-p2 ", clump_p2,
                 " --clump-r2 ", clump_r2,
                 " --clump-kb ", clump_kb,
                 " --out ", shQuote(fn, type = shell),
                 " --threads ", threads)
  system(fun2, ignore.stdout = T, ignore.stderr = T)
  res <- tryCatch({
    utils::read.table(paste(fn, ".clumped", sep = ""), header = T)
  }, error=function(e) {
    utils::read.table(paste(fn, sep = ""), header = T)
  }, warning=function(w) {
    utils::read.table(paste(fn, sep = ""), header = T)
  })
  unlink(paste(fn, "*", sep = ""))
  y <- subset(dat, !dat[["SNP"]] %in% res[["SNP"]])
  if (nrow(y) > 0) {
    message("Removing ", length(y[["SNP"]]), " of ", nrow(dat),
            " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["SNP"]] %in% res[["SNP"]]))
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
