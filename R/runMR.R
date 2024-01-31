## Func: run MR-IVW using TwoSampleMR based on TwoSampleMR::harmonise()
#' @title run MR-IVW using TwoSampleMR pkg
#'
#' @param datMR the harmonized data for MR analysis using TwoSampleMR::harmonise()
#'
#' @export
runMRIVW <- function(datMR) {
  datResIVW <- TwoSampleMR::mr(dat=datMR)
  datResIVW$pval <- format(datResIVW$pval,digits=3)
  return(datResIVW)
}


## Func: run MR-RAPS using mr.raps pkg based on TwoSampleMR::harmonise()
#' @title run MR-RAPS using TwoSampleMR pkg
#'
#' @param datMR the harmonized data for MR analysis using TwoSampleMR::harmonise()
#'
#' @export
runMRRAPS <- function(datMR) {
  datResult.MRRAPS <- with(
    datMR[datMR$mr_keep==TRUE,],
    mr.raps::mr.raps(b_exp=beta.exposure,
                     b_out = beta.outcome,
                     se_exp = se.exposure,
                     se_out = se.outcome))
  datRes <- data.frame(
    id.exposure = datMR$id.exposure[1],
    id.outcome = datMR$id.outcome[1],
    outcome = datMR$outcome[1],
    exposure = datMR$exposure[1],
    method = "MR RAPS",
    nsnp = sum(datMR$mr_keep == TRUE),
    b = datResult.MRRAPS$beta.hat,
    se = datResult.MRRAPS$beta.se,
    pval = format(datResult.MRRAPS$beta.p.value,digits=3)
  )
  return(datRes)
}


## Func: run MR-GRAPPLE using mr.divw pkg based on TwoSampleMR::harmonise()
#' @title run de-biased MR using TwoSampleMR pkg
#'
#' @param datMR the harmonized data for MR analysis using TwoSampleMR::harmonise()
#'
#' @export
runDebiasedMR <- function(datMR) {
  datDebiaedMR <- with(
    datMR[datMR$mr_keep == TRUE, ],
    mr.divw::mr.divw(beta.exposure = beta.exposure,
                     beta.outcome = beta.outcome,
                     se.exposure = se.exposure,
                     se.outcome = se.outcome,
                     lambda = 0))
  datRes <- data.frame(
    data.frame(
      id.exposure = datMR$id.exposure[1],
      id.outcome = datMR$id.outcome[1],
      outcome = datMR$outcome[1],
      exposure = datMR$exposure[1],
      method = "Debiased IVW",
      nsnp = sum(datMR$mr_keep == TRUE),
      b = datDebiaedMR$beta.hat,
      se = datDebiaedMR$beta.se,
      pval = format((1-stats::pnorm(abs(datDebiaedMR$beta.hat/datDebiaedMR$beta.se)))*2,
                    digits=3)))
  return(datRes)
}


## run MR-Clust using mrclust pkg based on TwoSampleMR::harmonise()
#' @title run MR-Clust analysis
#'
#' @param datMR the harmonized data for MR analysis using TwoSampleMR::harmonise()
#'
#' @export
#'
runMRClust <- function(datMR) {

  # select IVs used in MR
  datMRClust <- datMR[datMR$mr_keep==TRUE,]

  # calculate ratio-estimate and its standard error
  datMRClust$ratio.estimate <- with(datMRClust,beta.outcome/beta.exposure)
  datMRClust$se.ratio.estimate <- with(datMRClust,se.outcome/abs(beta.exposure))
  datMRClust$snpName <- with(datMRClust,
                             paste("chr",chr.exposure,":",pos.exposure,sep=""))
  # run an MR-Clust analysis
  datMRClust.Result <- with(
    datMRClust,
    mrclust::mr_clust_em(theta=ratio.estimate,
                         theta_se=se.ratio.estimate,
                         bx=beta.exposure,
                         by=beta.outcome,
                         bxse=se.exposure,
                         byse=se.outcome,
                         obs_name=SNP)
  )

  # Save scatter plot of the fitted model
  pltScatter <- datMRClust.Result$plots$two_stage +
    ggplot2::xlab(datMRClust$exposure[1]) +
    ggplot2::ylab(datMRClust$outcome[1])

  # select cluster with
  #     (1) at least 2 IVs
  #     (2) prob > 0.8
  datMRClust.best <- mrclust::pr_clust(
    datMRClust.Result$results$best,
    prob=0.8,min_obs=2)

  # update the scatter plot
  pltScatterProb08.2IVs <- with(
    datMRClust[datMRClust$SNP %in% datMRClust.best$observation,],
    mrclust::two_stage_plot(res=datMRClust.best,
                            bx=beta.exposure,
                            by=beta.outcome,
                            bxse=se.exposure,
                            byse=se.outcome,
                            obs_name=SNP) +
      ggplot2::xlab(datMRClust$exposure[1]) +
      ggplot2::ylab(datMRClust$outcome[1]))


  # search for traits associated with clusters of variants using PhenoScanner v2
  datMRClust.PhenoScanner <- with(
    datMRClust,
    tryCatch({
      phenoscanner::phenoscanner(snpquery=SNP,
                                 catalogue="GWAS",
                                 pvalue=5e-8,
                                 proxies="None",
                                 r2=0.8,
                                 build=37)
    }, error=function(e) {
      phenoscanner::phenoscanner(snpquery=SNP,
                                 catalogue="GWAS",
                                 pvalue=5e-8,
                                 proxies="None",
                                 r2=0.8,
                                 build=37)
    }))

  datMRClust.PheWAS <- datMRClust.PhenoScanner$results %>%
    dplyr::select(-rsid) %>%
    dplyr::left_join(datMRClust.PhenoScanner$snps[
      ,c("snp","afr","amr","eas","eur","sas","consequence",
         "protein_position","amino_acids","ensembl","hgnc")],
      by="snp") %>%
    dplyr::left_join(datMRClust.best[
      ,c("observation","cluster","probability")],
      by=c("snp"="observation"))

  return(list(datMRClust.all=datMRClust.Result$results$all,
              datMRClust.best=datMRClust.best,
              datMRClust.PheWAS=datMRClust.PheWAS,
              pltScatter=pltScatter,
              pltScatterProb08.2IVs=pltScatterProb08.2IVs))

}
