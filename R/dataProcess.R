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
#' @param clump_kb clumping window. By default, \code{clump_kb=1000}
#' @param clump_r2 clumping r2 cutoff. By default, \code{clump_r2=0.001}
#' @param threads number of threads used for plink
#'
#' @export
run_clump <- function(dat,
                      ldRef,
                      clump_kb = 1000,
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

#' @title Get LDmatrix for given SNPs concerning 1KG reference
#'
#' @param dat data frame from \code{TwoSampleMR::format_data()} for exposure GWAS
#' @param ldRef location of 1KG reference panel
#' @param pval pvalue threshold for selecting SNPs
#' @param threads number of threads used for plink
#' @param withAlleles indicator of whether LD matrix includes effect and other alleles
#'
#' @export
#'
getLDmatrix <- function(dat,ldRef,pval,threads=36,withAlleles=TRUE) {
  # ldMatrix <- ieugwasr::ld_matrix_local(
  #   variants = dat$SNP[dat$pval<pval],
  #   bfile = ldRef,
  #   plink_bin = plinkbinr::get_plink_exe()
  # )
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")
  fn <- tempfile()
  utils::write.table(data.frame(SNP=dat[["SNP"]],P=dat[["pval.exposure"]]),
                     file = fn, row.names = F, col.names = T, quote = F)
  fun1 <- paste0(shQuote(plinkbinr::get_plink_exe(), type = shell),
                 " --bfile ", shQuote(ldRef, type = shell),
                 " --extract ", shQuote(fn, type = shell),
                 " --make-just-bim ",
                 " --keep-allele-order ",
                 " --out ", shQuote(fn, type = shell),
                 " --threads ", threads)
  system(fun1)
  bim <- utils::read.table(paste0(fn, ".bim"), stringsAsFactors = FALSE)
  fun2 <- paste0(shQuote(plinkbinr::get_plink_exe(), type = shell),
                 " --bfile ", shQuote(bfile, type = shell),
                 " --extract ", shQuote(fn, type = shell),
                 " --r square ",
                 " --keep-allele-order ",
                 " --out ", shQuote(fn, type = shell),
                 " --threads ", threads)
  system(fun2)
  ldMatrix <- utils::read.table(paste0(fn, ".ld"), header = FALSE) %>% as.matrix()
  snp_list <- bim$V2
  if (withAlleles) {
    rownames(ldMatrix)<-colnames(ldMatrix)<-paste(bim$V2,bim$V5,bim$V6,sep="_")
  } else {
    rownames(ldMatrix) <- colnames(ldMatrix) <- bim$V2
  }
  # snp_list <- unlist(lapply(
  #   row.names(ldMatrix),
  #   FUN=function(i) {
  #     strsplit(i,"_")[[1]][1]
  #   }))
  dat <- dat[dat$SNP %in% snp_list,]
  idx <- unlist(lapply(snp_list, FUN=function(i) which(dat$SNP %in% i)))
  dat <- dat[idx,]
  return(list(LDmatrix=ldMatrix,dat=dat))
}


run_SuSiE <- function(dat,ldMatrix,nSampleSize,xlabel) {
  ## calculate the z-score for fine-mapping
  dat$zscore <- with(datGeneMR,beta.exposure/se.exposure)
  fittedSuSiE <- with(
    datGeneMR,
    susieR::susie_rss(z=zscore,
                      R=ldMatrix,
                      L=10,n=nSampleSize,
                      estimate_residual_variance=TRUE,
                      estimate_prior_variance=TRUE))
  calPIP <- summary(fittedSuSiE)$vars
  selectedSNPs <- sort(calPIP$variable[calPIP$cs>0])

  # susieR::susie_plot(fittedSuSiE, y="PIP",
  #                    b=datGeneStat$beta.exposure[1:00],
  #                    xlab=paste(uniprot_gn_symbol,
  #                               " (Chr",
  #                               datCodingGene$chrpos[datCodingGene$uniprot_gn_symbol==uniprot_gn_symbol][1],
  #                               ")",sep=""))
  #
  # gaston::LD.plot(ldSNP,snp.positions=datGeneStat$pos.exposure,
  #                 polygon.par = list(border = NA),
  #                 write.ld = NULL)

  return(list(PIPs=calPIP,credibleSet=selectedSNPs, ## fine-mapping outputs
              datMR=datMR,ldSNP=ldSNP,              ## data set for MR
              datMRResult=datMRResult))             ## MR results with correlated IVs
}


#' @title remove ".exposure" and ".outcome" from a data frame
#'
#' @param dat GWAS summary statistics from \code{TwoSampleMR::format_data()}
#'
#' @export
rnDat <- function(dat) {
  names(dat) <- gsub(".exposure","",names(dat))
  names(dat) <- gsub(".outcome","",names(dat))
  return(dat)
}


#' @title calculate number of the overlapped SNPs between two traits given a pval threshold
#'
#' @param dat1 GWAS summary statistics from \code{TwoSampleMR::format_data()} for trait 1
#' @param dat2 GWAS summary statistics from \code{TwoSampleMR::format_data()} for trait 2
#' @param pval pvalue threshod for selecting the overlapped SNPs
#'
#' @export
checkSNPsOverlap <- function(dat1,dat2,pval) {
  dat1 <- rnDat(dat1)
  dat2 <- rnDat(dat2)
  snps <- intersect(dat1$SNP[dat1$pval<pval],dat2$SNP[dat2$pval<pval])
  return(list(nSNPs=length(snps),SNPs=snps))
}

#' @title get the overlapped SNPs between two traits for colocalization analyses
#'
#' @param dat1 GWAS summary statistics from \code{TwoSampleMR::format_data()} for trait1
#' @param trait1 trait name for dat1
#' @param dat2 GWAS summary statistics from \code{TwoSAmpleMR::format_data()} for trait2
#' @param trait2 trait name for dat2
#' @param ldRef 1KG reference panel in plink format
#' @param pval pvalue threshold for selecting SNPs with \code{pvalue<pval}
#'
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
