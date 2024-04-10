# drugTargetScreen::dat4UKBNealeLab() is more efficient than this function(
# 
# this.path::this.dir()

library(data.table)
library(ggplot2)
library(ggpubr)

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


## plt scatter plot of MR
plt.Scatter <- function(datMR,title,legend.idx=TRUE) {
    if(legend.idx) {
  TwoSampleMR::mr_scatter_plot(TwoSampleMR::mr(datMR),dat=datMR)[[1]] +
         ggplot2::labs(x="Genetic effects on LDL-C",
                       y="Genetic effects on CAD",
                       title=title) +
         ggplot2::theme_classic() +
         ggplot2::theme(legend.position=c(0.2,0.8),
                        legend.justification="left",
                        legend.box.just="top",
                        # legend.direction="vertical",
                        axis.line=element_line(linetype = "solid"),
			            panel.background=element_blank(),
			            panel.border=element_blank(),
			            panel.grid.major=element_blank(),
			            panel.grid.minor=element_blank(),
			            plot.background=element_blank(),
			            plot.title=element_text(size=11,hjust=0.5)) +
         ggplot2::guides(color=ggplot2::guide_legend(ncol=2,byrow=TRUE)) +
         ggsci::scale_color_jama()
    } else {
          TwoSampleMR::mr_scatter_plot(TwoSampleMR::mr(datMR),dat=datMR)[[1]] +
         ggplot2::labs(x="Genetic effects on LDL-C",
                       y="Genetic effects on CAD",
                       title=title) +
         ggplot2::theme_classic() +
         ggplot2::theme(legend.position="none",
                        # legend.justification="right",
                        # lgend.box.just="bottom",
                        # legend.direction="vertical",
                        axis.line=element_line(linetype = "solid"),
			            panel.background=element_blank(),
			            panel.border=element_blank(),
			            panel.grid.major=element_blank(),
			            panel.grid.minor=element_blank(),
			            plot.background=element_blank(),
			            plot.title=element_text(size=11,hjust=0.5)) +
         ggplot2::guides(color=ggplot2::guide_legend(ncol=2,byrow=TRUE)) +
         ggsci::scale_color_jama()
    }

}

## plt volcano plot
plt.Volcano <- function(mainMethod,gene) {
  ggplot(data=MRResultPheMR[method %in%  mainMethod & external_gene_name %in% gene],
                      aes(x=b,y=-log10(as.numeric(pval)))) +
    geom_jitter(aes(colour= factor(dir)),show.legend=FALSE) +
    ggrepel::geom_text_repel(
      data=MRResultPheMR[idx==TRUE & method %in% mainMethod & external_gene_name %in% c("CETP")],
      aes(label=sidePhe),box.padding=0.5,max.overlaps=Inf,
      min.segment.length=0) +
    ggsci::scale_color_jama() + ggsci::scale_fill_jama() +
    ggplot2::geom_vline(aes(xintercept=log(1.2)),linetype=2,color="gray") +
    ggplot2::geom_vline(aes(xintercept=log(0.8)),linetype=2,color="gray") +
    geom_hline(ggplot2::aes(yintercept=-log10(FDR.level/78)),color="gray",linetype=2) +
    scale_x_continuous(name="log(OR)",limits=c(-4,4)) +
    scale_y_continuous(name=expression("-log"[10]*"(p-value)")) +
    theme_classic() +
    theme(legend.position=NULL,
          axis.line = element_line(linetype = "solid"),
			    panel.background=element_blank(),
			    panel.border=element_blank(),
			    panel.grid.major=element_blank(),
			    panel.grid.minor=element_blank(),
			    plot.background=element_blank(),
			    plot.title=element_text(size=11,hjust=0.5)) +
    ggtitle(paste("Potential side effects of drug-target on ",gene,sep=""))
}


## plt CIF using object obtained from gform
plt.CIF <- function(gform.fit,title,legend.idx=TRUE){
    dat <- melt(
        gform.fit$result,
        id.vars=c("k","Interv."),
        measure.vars=c("NP Risk","g-form risk"),
        variable.name="Method",
        value.name="CIF")
    if (legend.idx) {
    pltCIF <- ggplot(dat,aes(x=k,y=CIF,col=Method)) +
        geom_line() + #geom_point() +
        ggsci::scale_color_jama() + ggsci::scale_fill_jama() +
        theme_classic() +
        theme(legend.position=c(0.2,0.8),
              legend.justification="top",
              legend.box.just="left",
              # legend.direction="vertical",
              legend.title=element_blank(),
              axis.line = element_line(linetype = "solid"),
			        panel.background=element_blank(),
			        panel.border=element_blank(),
			        panel.grid.major=element_blank(),
			        panel.grid.minor=element_blank(),
			        plot.background=element_blank(),
			        plot.title=element_text(size=11,hjust=0.5)) +
        scale_x_continuous(name="Follow-up time (yrs)", limits=c(0,30)) +
        scale_y_continuous(name="Cumulative risk",limits=c(0,0.3)) +
        labs(title=title)
    } else {
    pltCIF <- ggplot(dat,aes(x=k,y=CIF,col=Method)) +
        geom_line() + #geom_point() +
        ggsci::scale_color_jama() + ggsci::scale_fill_jama() +
        theme_classic() +
        theme(legend.position="none",
              axis.line = element_line(linetype = "solid"),
			        panel.background=element_blank(),
			        panel.border=element_blank(),
			        panel.grid.major=element_blank(),
			        panel.grid.minor=element_blank(),
			        plot.background=element_blank(),
			        plot.title=element_text(size=11,hjust=0.5)) +
        scale_x_continuous(name="Follow-up time (yrs)", limits=c(0,30)) +
        scale_y_continuous(name="Cumulative risk",limits=c(0,0.3)) +
        labs(title=title)      
    }
    return(pltCIF)
}
## plt CIF for hypothetical interventions
plt.Int <- function(gform.fit=gform_ASCVD,title,int.lab,legend.idx=TRUE) {
    dat <- melt(
        gform.fit$result,
        id.vars=c("k","Interv."),
        measure.vars=c("g-form risk"),
        variable.name="Method",
        value.name="CIF")
    dat <- dat[!is.na(dat$CIF) & dat$`Interv.` %in% c(0,1,2),]
    dat$Int <- factor(dat$`Interv.`,levels=0:2,labels=int.lab[1:3])
    if(legend.idx) {
    pltInt <- ggplot(dat,aes(x=k,y=CIF,col=Int)) +
        geom_line() + #geom_point() +
        ggsci::scale_color_jama() + ggsci::scale_fill_jama() +
        theme_classic() +
        theme(legend.position=c(0.2,0.8),
              legend.justification="top",
              legend.box.just="left",
              # legend.direction="vertical",
              legend.title=element_blank(),
              axis.line = element_line(linetype = "solid"),
			        panel.background=element_blank(),
			        panel.border=element_blank(),
			        panel.grid.major=element_blank(),
			        panel.grid.minor=element_blank(),
			        plot.background=element_blank(),
			        plot.title=element_text(size=11,hjust=0.5)) +
        scale_x_continuous(name="Follow-up time (yrs)", limits=c(0,30)) +
        scale_y_continuous(name="Cumulative risk",limits=c(0,0.3)) +
        labs(title=title)
    } else {
      pltInt <- ggplot(dat,aes(x=k,y=CIF,col=Int)) +
        geom_line() + #geom_point() +
        ggsci::scale_color_jama() + ggsci::scale_fill_jama() +
        theme_classic() +
        theme(legend.position="none",
              axis.line = element_line(linetype = "solid"),
			        panel.background=element_blank(),
			        panel.border=element_blank(),
			        panel.grid.major=element_blank(),
			        panel.grid.minor=element_blank(),
			        plot.background=element_blank(),
			        plot.title=element_text(size=11,hjust=0.5)) +
        scale_x_continuous(name="Follow-up time (yrs)", limits=c(0,30)) +
        scale_y_continuous(name="Cumulative risk",limits=c(0,0.3)) +
        labs(title=title)
    }
    return(pltInt)
}
format95CI <- function(x,lci,uci) {
    paste(sprintf("%.2f",x), " (",
          sprintf("%.2f",lci), " to ",
          sprintf("%.2f",uci), ")", sep="")
}
getRiskRatio <- function(gform.fit=gform_ASCVD,int.lab) {
    dat <- gform.fit$result
    dat$Int <- factor(dat$`Interv.`,levels=0:3,labels=int.lab)
    dat <- dat[,.SD[.N],by=.(Int)]
    dat <- dat[,.(`NP Risk`=sprintf("%.3f",`NP Risk`),
                  `g-form risk`=sprintf("%.3f",`g-form risk`),
                  RiskRatio=format95CI(`Risk ratio`,`RR lower 95% CI`,`RR upper 95% CI`),
                  RiskDiff=format95CI(`Risk difference`,`RD lower 95% CI`,`RD upper 95% CI`))
               ,by=.(Int)]
}

getMRestimates <- function(MRResult,mainMethod) {
  data.table::setDT(MRResult)
  MRResult <- MRResult[
    ,`:=`(OR=format95CI(exp(b),exp(b-1.96*se),exp(b+1.96*se)))
    ,by=.(method)]
  data.table::setDF(MRResult)
  return(MRResult[MRResult$method %in% mainMethod,])  
}


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