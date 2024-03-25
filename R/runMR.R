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
#' @importFrom magrittr %>%
#'
#' @export
runMRClust <- function(datMR) {

  rsid <- NULL

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


  # # search for traits associated with clusters of variants using PhenoScanner v2
  # datMRClust.PhenoScanner <- with(
  #   datMRClust,
  #   tryCatch({
  #     phenoscanner::phenoscanner(snpquery=SNP,
  #                                catalogue="GWAS",
  #                                pvalue=5e-8,
  #                                proxies="None",
  #                                r2=0.8,
  #                                build=37)
  #   }, error=function(e) {
  #     phenoscanner::phenoscanner(snpquery=SNP,
  #                                catalogue="GWAS",
  #                                pvalue=5e-8,
  #                                proxies="None",
  #                                r2=0.8,
  #                                build=37)
  #   }))
  #
  # datMRClust.PheWAS <- datMRClust.PhenoScanner$results %>%
  #   dplyr::select(-rsid) %>%
  #   dplyr::left_join(datMRClust.PhenoScanner$snps[
  #     ,c("snp","afr","amr","eas","eur","sas","consequence",
  #        "protein_position","amino_acids","ensembl","hgnc")],
  #     by="snp") %>%
  #   dplyr::left_join(datMRClust.best[
  #     ,c("observation","cluster","probability")],
  #     by=c("snp"="observation"))

  return(list(datMRClust.all=datMRClust.Result$results$all,
              datMRClust.best=datMRClust.best,
              # datMRClust.PheWAS=datMRClust.PheWAS,
              pltScatter=pltScatter,
              pltScatterProb08.2IVs=pltScatterProb08.2IVs))

}


## run MR-IVW, MR-RAPS and debiased-MR in all based on TwoSampleMR::harmonise()
#' @title run MR analyses based on MR-Clust
#'
#' @param datMR data frame from TwoSampleMR::harmonised_data()
#' @param datMRClust.best data frame obtained from MR-Clust with best selection
#'
#' @export
runMRClustBest <- function(datMR,datMRClust.best) {
  lClust <- unique(datMRClust.best$cluster)
  datResult <- do.call(
    "rbind",
    lapply(lClust, function(i) {
      dat <- datMR[
        datMR$SNP %in% datMRClust.best$observation[datMRClust.best$cluster==i],]
      datRes <- rbind(runMRIVW(datMR=dat),
                      runMRRAPS(datMR=dat),
                      runDebiasedMR(datMR=dat))
      datRes$MRCluster <- paste("Cluster ",i,sep="")
      datRes
    }))
  return(datResult[,!names(datResult) %in% c("id.exposure","id.outcome")])
}


#' @title run drug target (cis-) MR analyses with annotated data
#'
#' @param datExp data frame for annotated exposure GWAS using \code{annotateGeneName()}
#' @param datMR data frame from \code{TwoSampleMR::harmonised_data()}
#' @param FDR.level threshold for controlling false positive rate using q-value
#'
#' @export
runCisMR <- function(datExp,datMR,FDR.level=0.01) {
  pval <- b <- GENCODE_name <- NULL
  uniCodingGene <- unique(datExp$uniprot_gn_symbol)
  uniCodingGene <- uniCodingGene[!uniCodingGene %in% c(NA,"")]
  datCisMRResult <- do.call(
    "rbind",
    lapply(uniCodingGene, function(i) {
      cat(paste("Coding GENE: ",i," | ",sep=""))
      datCisMR <- datMR[datMR$SNP %in% datExp$SNP[datExp$uniprot_gn_symbol %in% i],]
      if (nrow(datCisMR) > 0) {
        tryCatch({
          datRes <- runMRIVW(datMR = datCisMR)
          datRes$GENCODE_name <- i
          datRes
        }, error=function(e) NULL)
      }
    }))
  datCisMRResult <- datCisMRResult[,!names(datCisMRResult) %in% c("id.exposure","id.outcome")]
  ## Clear results
  datCisMRResult <- within(datCisMRResult,{
    qval = qvalue::qvalue(p=as.numeric(pval),fdr.level=FDR.level)$qvalues
    idx = qval < FDR.level
    dir = ifelse(b<0,"Protective","Harzard")
  })
  ## volcano plot
  pltVolcano <- ggplot2::ggplot(data=datCisMRResult,
                       ggplot2::aes(x=b,y=-log10(as.numeric(pval)))) +
    ggplot2::geom_jitter(ggplot2::aes(colour= factor(dir)),
                show.legend=FALSE) +
    ggrepel::geom_text_repel(
      data=datCisMRResult[datCisMRResult$idx==TRUE,],
      ggplot2::aes(label=GENCODE_name),box.padding=0.5,max.overlaps=Inf,
      min.segment.length=0) +
    ggsci::scale_color_jama() + ggsci::scale_fill_jama() +
    ggplot2::geom_hline(ggplot2::aes(yintercept=-log10(FDR.level)),
               color="gray",linetype=2) +
    ggplot2::scale_x_continuous(name="log(OR)") +
    ggplot2::scale_y_continuous(name=expression("-log"[10]*"(p-value)")) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position=NULL) +
    ggplot2::ggtitle(paste("Cis-MR with q-value < ",FDR.level,sep=""))

  return(list(datCisMRResult=datCisMRResult,pltVolcano=pltVolcano))
}


#' @title run Gene-based MR analyses with correlated IVs
#'
#' @param datMR data frame from \code{TwoSampleMR::harmonised_data()}
#' @param ldRef 1KG reference panel for getting LD matrices
#' @param pval threshold for selecting IVs at p-value
#' @param corr \code{corr=c(TRUE,FALSE)} indicator of inclusion of correlated IVs or not
#'
#' @export
runGeneMendelianRandomization <- function(datMR,ldRef,pval=1e-4,corr=TRUE) {
  SNP <- pval.exposure <- NULL
  ## Pipeline:
  ## (1) calculate the LD matrix for the specific gene region
  ## (2) clean the selected IVs for the specific gene region
  ## (3) prepare the data used for MendelianRandomization pkg
  ## (4) run MR analyses, including
  ##      (4.0) penalized inverse-variance weighted method
  ##      (4.1) debiased MR
  ##      (4.2) mr-median
  ##      (4.3) mr-mode
  ##      (4.4) mr-pcgmm: particularly useful for gene-level MR
  ## Calculate the LD matrix for MR analysis with correlated IVs
  if (corr) {
    ldMatrix <- ieugwasr::ld_matrix_local(
      variants = datMR$SNP[datMR$mr_keep & datMR$pval.exposure<pval],
      bfile = ldRef,
      plink_bin = plinkbinr::get_plink_exe()
    )
    snp_list <- unlist(lapply(
      row.names(ldMatrix),
      FUN=function(i) {
        strsplit(i,"_")[[1]][1]
      }))
    datMR <- datMR[datMR$SNP %in% snp_list,]
    idx <- unlist(lapply(snp_list, FUN=function(i) which(datMR$SNP %in% i)))
    datMR <- datMR[idx,]
    datMRNew <- with(
      datMR,
      MendelianRandomization::mr_input(
        bx=beta.exposure,
        bxse=se.exposure,
        by=beta.outcome,
        byse=se.outcome,
        correlation=ldMatrix))
  } else {
    datMRNew <- with(
      datMR,
      MendelianRandomization::mr_input(
        bx=beta.exposure,
        bxse=se.exposure,
        by=beta.outcome,
        byse=se.outcome))
  }

  ## run Mendelian randomization
  return(runMendelianRandomization(datMRInput=datMRNew))
}


#' @title run MR using \pkg{MendelianRandomization}
#'
#' @param datMRInput data frame from \code{MendelianRandomization::mr_input()}
#'
#' @export
runMendelianRandomization <- function(datMRInput) {
  ## run Mendelian randomization
  datMR_pivw <- attributes(MendelianRandomization::mr_pivw(datMRInput))
  datMR_divw <- attributes(MendelianRandomization::mr_divw(datMRInput))
  datMR_median<-attributes(MendelianRandomization::mr_median(datMRInput))
  datMR_mode <- attributes(MendelianRandomization::mr_mbe(datMRInput))
  # datMR_pcgmm<- attributes(MendelianRandomization::mr_pcgmm(datMRInput,nx=nx,ny=ny))

  datMRResult <- rbind(
    data.frame(method="Penalised IVW",b=datMR_pivw$Estimate,se=datMR_pivw$StdError,cond=datMR_pivw$Condition),
    data.frame(method="Debiased IVW",b=datMR_divw$Estimate,se=datMR_divw$StdError,cond=datMR_divw$Condition),
    data.frame(method="Weighted median",b=datMR_median$Estimate,se=datMR_median$StdError,cond=NA),
    data.frame(method="Simple Mode",b=datMR_mode$Estimate,se=datMR_mode$StdError,cond=NA)#,
    # data.frame(method="PCA of GMM",b=datMR_pcgmm$Estimate,se=datMR_pcgmm$StdError,cond=NA)
  )
  datMRResult$pval <- with(datMRResult,
                           format(2*pnorm(abs(b/se),lower.tail=FALSE),digits=3))
  return(datMRResult)
}



#' @title run MR analysis based on SuSiER
#'
#' @param datGeneMR data frame with gene-level data from \code{TwoSampleMR::harmonised_data()}
#' @param ldRef 1KG reference panel for getting LD matrices
#' @param clump_r2 correlation in LD. By default, \code{clump_r2=0.2}
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr .
#'
#' @export
runSuSiERMR <- function(datGeneMR,ldRef,clump_r2=0.2) {
  SNP <- pval.exposure <- NULL
  ## clumping
  retain_snps <- datGeneMR %>%
    dplyr::select(rsid=SNP,pval=pval.exposure) %>%
    ieugwasr::ld_clump_local(
      .,
      bfile = ldRef,
      clump_kb = 1000,
      clump_r2 = clump_r2,
      clump_p = 5e-8,
      plink_bin = plinkbinr::get_plink_exe()) %>%
    {.$rsid}
  ## get the ld matrix
  ldMatrix <- ieugwasr::ld_matrix_local(
    variants = datGeneMR$SNP,
    bfile = ldRef,
    plink_bin = plinkbinr::get_plink_exe()
  )
  snp_list <- unlist(lapply(
    row.names(ldMatrix),
    FUN=function(i) {
      strsplit(i,"_")[[1]][1]
    }))
  datGeneMR <- datGeneMR[datGeneMR$SNP %in% snp_list,]
  idx <- unlist(lapply(snp_list, FUN=function(i) which(datGeneMR$SNP %in% i)))
  datGeneMR <- datGeneMR[idx,]
  ## calculate the z-score for fine-mapping
  datGeneMR$zscore <- with(datGeneMR,beta.exposure/se.exposure)
  fittedSuSiE <- with(
    datGeneMR,
    susieR::susie_rss(z=zscore,
                      R=ldMatrix,
                      L=10,n=10000,
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

  ## run fine-mapping based MR
  datMR <- datGeneMR[selectedSNPs,]
  ldSNP <- ldMatrix[selectedSNPs,selectedSNPs]
  datMRNew <- with(
    datMR,
    MendelianRandomization::mr_input(
      bx=beta.exposure,
      bxse=se.exposure,
      by=beta.outcome,
      byse=se.outcome,
      correlation=ldSNP))

  datMRResult <- runMendelianRandomization(datMRInput=datMRNew)
  return(list(PIPs=calPIP,credibleSet=selectedSNPs, ## fine-mapping outputs
              datMR=datMR,ldSNP=ldSNP,              ## data set for MR
              datMRResult=datMRResult))             ## MR results with correlated IVs

}



