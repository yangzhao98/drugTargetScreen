# devtools::install_github("yangzhao98/drugTargetScreen",force=TRUE)
## setup path info
path.C <- "D:"
# path.C <- "D:"
path.ensembl <- paste0(path.C,"/PublicData/Ensembl")
ldRef <- paste0(path.C,"/PublicData/1KG_v3_OpenGWAS/EUR")
path.gclc <- paste0(path.C,"/PublicData/GLGC/GLGClipids2021")
path.ieuGWAS <- paste0(path.C,"/PublicData/ieuGWAS")
path.cad <- paste0(path.C,"/PublicData/CARDIoGRAMplusC4D")
path.ukb <- paste0(path.C,"/PublicData/UKB_nealelab")
path.PheMR <- paste0(path.C,"/PublicData/UKB_nealelab")
path.gtexv8 <- paste0(path.C,"/PublicData/GTEx_v8/eQTLGen_phase1") 
path.gtexv8.organ <- paste0(path.C,"/PublicData/GTEx_v8/GTEx_Analysis_v8_eQTL")

## load GRCh37 with gene info
## ---
## GRCh37 gene info
posCtrl <- c("APOB","LDLR","LPA")
datGRCh37 <- data.table::fread(
  paste0(path.ensembl,"/",
        grep("GRCh37.87.geneName",list.files(path.ensembl),value="TRUE")))
datPosCtrlNew <- datGRCh37[external_gene_name %in% posCtrl][
  ,.SD[.N],by=.(external_gene_name)][order(chromosome_name)][
  external_gene_name != "CETP"
]
datPosCtrlNew
