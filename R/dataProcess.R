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


#' @title clumping the data
#'
#' @param dat data frame with at least two columns names \code{c("SNP","pval.exposure")}
#' @param plink_exe location of plink.exe. By default, \code{plink_exe=plinkbinr::get_plink_exe()}
#' @param ldRef 1KG reference panel with .bed format
#' @param clump_p1 p-value threshold for the index SNP. By default, \code{clump_p1=1}
#' @param clump_p2 p-value threshold for the secondary SNP. By default, \code{clump_p2=1}
#' @param clump_kb clumping window. By default, \code{clump_kb=10000}
#' @param clump_r2 clumping r2 cutoff. By default, \code{clump_r2=0.001}
#' @param tempdir temporal direct for storing plink data set
#'
#' @export
plink_clump <- function(dat,
                        plink_exe,
                        ldRef,
                        clump_kb = 10000,
                        clump_r2 = 0.001,
                        clump_p1 = 1,
                        clump_p2 = 1,
                        tempdir = "temp") {
  # Make textfile
  snps <- dat$SNP
  if (!("pval" %in% colnames(dat))) dat$pval <- dat$pval.exposure
  pvals <- pmin(1, dat$pval)
  dir.create(file.path(tempdir), showWarnings = FALSE)

  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  fn <- tempfile(tmpdir = tempdir)
  require(data.table)
  data.table::fwrite(data.frame(SNP=snps, P=pvals), file=fn, row=F, col=T, qu=F, sep = " ")

  fun2 <- paste0(
    # shQuote(plink_exe),
    plink_exe,
    #  shQuote("plink"),
    #    " --bfile ", shQuote(refdat, type=shell),
    #    " --clump ", shQuote(fn, type=shell),
    " --bfile ", ldRef,
    " --clump ", fn,
    " --clump-p1 ", clump_p1,
    " --clump-p2 ", clump_p2,
    " --clump-r2 ", clump_r2,
    " --clump-kb ", clump_kb,
    " --out ", shQuote(fn, type=shell)
  )
  system(fun2, ignore.stdout = T, ignore.stderr = T)
  # a <- fread(paste(fn, ".clumped", sep=""), he=T)
  a <- data.table::fread(paste(fn, sep=""), he=T)
  unlink(paste(fn, "*", sep=""))
  a <- a[, c(3, 5)]
  a$temp.p <- round(log10(a$P))
  dat$temp.p <- round(log10(dat$pval))
  a <- merge(a, dat, by = c("SNP", "temp.p"))

  a <- a[, -(2:3)]

  unlink(tempdir, recursive = T)

  return(a)
}
