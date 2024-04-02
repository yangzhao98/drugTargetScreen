#' @title Calculate the area under a time-series curve
#' @param t the time
#' @param x the value at each time
#' @export
area <- function(t,x) sum(diff(t)*(x[-length(x)]+x[-1])*0.5)

#' @title Import data files with SPSS format
#' @param name The name of the SPSS data files
#' @param path The path of the SPSS data files
#' @importFrom foreign read.spss
#' @export
loadSPSS <- function(name,path) {
  data.table::as.data.table(foreign::read.spss(
    paste(path,name,sep=""),
    trim_values=TRUE,
    stringsAsFactors=FALSE,
    add.undeclared.levels="no"))
}

