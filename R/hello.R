# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

# library(devtools)
# library(roxygen2)
# library(this.path)
#' @title get system-specific paths
#'
#' @param path a specific path without root fold info
#'
#' @export
getPaths <- function(path) {
  if(!Sys.info()["sysname"] %in% "Windows") {
    paste("/mnt/c",path,sep="")
  } else {
    paste("C:",path,sep="")
  }
}
