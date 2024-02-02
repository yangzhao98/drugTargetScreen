#' @title read harmonized summary statistics from GWAS Catalog
#'
#' @param filename location and file name of a specific data set
#' @param type type of the data set, \code{type=c("exposure","outcome")}
#'
#' @export

readHarmonizedGWASCatalogGRCh37 <- function(filename,type,...) {
  type <- match.arg(type,c("exposure","outcome"))
  TwoSampleMR::format_data(
    dat=data.table::setDF(data.table::fread(filename,...)),
    type=type,
    snp_col="hm_rsid",
    beta_col="hm_beta",
    se_col="standard_error",
    eaf_col="hm_effect_allele_frequency",
    effect_allele_col="hm_effect_allele",
    other_allele_col="hm_other_allele",
    pval_col="p_value",
    chr_col="hm_chrom",
    pos_col="hm_pos")
}


# datExp <- with(datExp, datExp[!is.na(SNP) & pval.exposure<1e-4,])

#' @title read summary statistics from deCODE using Olink Explore 1536 Panel
#'
#' @param filename location and file name of a specific data set
#' @param type type of the data set, \code{type=c("exposure","outcome")}
#'
#' @export
readdeCODEProteomicsGRCh38 <- function(filename,type) {
  type <- match.arg(type,c("exposure","outcome"))
  TwoSampleMR::format_data(
    dat=data.table::setDF(data.table::fread(filename)),
    type=type,
    snp_col="rsids",
    beta_col="Beta",
    se_col="SE",
    eaf_col="ImpFreqA1",
    effect_allele_col="A1",
    other_allele_col="A0",
    pval_col="Pval",
    chr_col="Chrom",
    pos_col="Pos",
    samplesize_col="N")
}
