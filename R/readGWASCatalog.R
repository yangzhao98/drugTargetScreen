## Assign the phyiscal path for GWASCatalog
immuneTraitsPath <- "C:/PublicData/NatGenet2020_ImmunityTraits/harmonizedSummaryStats"
sys_info <- Sys.info()["sysname"]
if (!sys_info %in% "Windows") immuneTraitsPath <- gsub("C:","/mnt/c",immuneTraitsPath)

## Get list of all data files
datFiles <- list.files(immuneTraitsPath)
datFiles <- grep(".tsv.gz",datFiles,value=TRUE)

# datTmp <- data.table::fread(paste(immuneTraitsPath,"/",datFiles[1],sep=""))

## Annotation with rsID and geneInfo
# Get GRCh37 info
GRCh37RefPath <- "C:/PublicData/Ensembl/Ensembl2Position.GRCh37.87.geneName.txt"

#
datExp <- TwoSampleMR::format_data(
  dat=data.table::setDF(
    data.table::fread(paste(immuneTraitsPath,"/",datFiles[1],sep=""))),
  type="exposure",
  snp_col="variant_id",
  beta_col="hm_beta",
  se_col="standard_error",
  eaf_col="hm_effect_allele_frequency",
  effect_allele_col="hm_effect_allele",
  other_allele_col="hm_other_allele",
  pval_col="p_value",
  chr_col="hm_chrom",
  pos_col="hm_pos",
  samplesize_col="n")
