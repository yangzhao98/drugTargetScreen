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
  if(sum(names(dat) %in% c("pval.exposure"))==1) {
    utils::write.table(data.frame(SNP=dat[["SNP"]],P=dat[["pval.exposure"]]),
                       file = fn, row.names = F, col.names = T, quote = F)
  } else if (sum(names(dat) %in% c("pval.outcome"))==1) {
    utils::write.table(data.frame(SNP=dat[["SNP"]],P=dat[["pval.outcome"]]),
                       file = fn, row.names = F, col.names = T, quote = F)
  } else {
    utils::write.table(data.frame(SNP=dat[["SNP"]],P=dat[["pval"]]),
                       file = fn, row.names = F, col.names = T, quote = F)
  }

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
                 " --bfile ", shQuote(ldRef, type = shell),
                 " --extract ", shQuote(fn, type = shell),
                 " --r square ",
                 " --keep-allele-order ",
                 " --out ", shQuote(fn, type = shell),
                 " --threads ", threads)
  system(fun2)
  ldMatrix <- as.matrix(utils::read.table(paste0(fn, ".ld"), header = FALSE))
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

