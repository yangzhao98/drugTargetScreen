.onLoad <- function(libname=find.package("drugTargetScreen"),pkgname="drugTargetScreen"){

  # CRAN Note avoidance
  if (getRversion() >= "3.1.0")
    utils::globalVariables(c("%>%", "."))
  invisible()
}
