.onLoad <- function(libname=find.package("drugTargetScreen"),pkgname="drugTargetScreen"){

  # CRAN Note avoidance
  if (getRversion() >= "3.1.0")
    utils::globalVariables(c("%>%", ".", "bfile", "plink_bin","pdf",":=",
                             "n_complete_samples","ra","se","ea","variant",
                             "expected_case_minor_AC","minor_AF",
                             "minor_allele","nCase"))
  invisible()

}
