#' @title plot Manhattan plot
#'
#' @param datManhattan data frame with the first three column being \code{c("SNP","Chromosome","Position")} with the following columns being \code{pval} for traits of interest
#' @param leadingSNPs data frame with at least \code{c("SNP","uniprot_gn_symbol")} info
#' @param filename location and name of the figure
#' @param type the type of the Manhattan plot, with \code{type=c("circle","manhattan")}
#'
#' @export
pltManhattan <- function(datManhattan,leadingSNPs,filename,
                         type=c("circle","manhattan")) {
  type <- match.arg("type",c("circle","manhattan"))
  CMplot::CMplot(datManhattan,
                 type="p",plot.type=type,LOG10=TRUE,
                 highlight=list(leadingSNPs$SNP,NULL,NULL),
                 highlight.text=list(leadingSNPs$uniprot_gn_symbol,NULL,NULL),
                 # highlight.type="l",
                 highlight.cex=0.5,
                 threshold=5e-8,threshold.col="black",
                 col=c("gray40","gray80"),
                 cex=0.8,
                 file="png",file.name=filename,
                 file.output=TRUE,multracks=TRUE,verbose=TRUE,dpi=300)
}

