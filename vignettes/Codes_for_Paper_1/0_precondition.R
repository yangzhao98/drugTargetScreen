# devtools::install_github("yangzhao98/drugTargetScreen",force=TRUE)

## setup path info
path.C <- "C:"
path.ensembl <- paste(path.C,"/PublicData/Ensembl",sep="")
ldRef <- paste(path.C,"/PublicData/1KG_v3_OpenGWAS/EUR",sep="")
path.gclc <- paste(path.C,"/PublicData/GLGC/GLGClipids2021",sep="")
path.cad <- paste(path.C,"/PublicData/CARDIoGRAMplusC4D",sep="")
path.ukb <- paste(path.C,"/PublicData/UKB_nealelab",sep="")
path.PheMR <- paste(path.C,"/PublicData/UKB_nealelab",sep="")
path.gtexv8 <- paste(path.C,"/PublicData/GTEx_v8/eQTLGen_phase1",sep="") 

## load GRCh37 with gene info
## ---
## GRCh37 gene info
datGRCh37 <- data.table::fread(
  paste(path.ensembl,"/",
        grep("GRCh37.87.geneName",list.files(path.ensembl),value="TRUE"),
        sep="")
)
datPosCtrlNew <- datGRCh37[external_gene_name %in% posCtrl][
  ,.SD[.N],by=.(external_gene_name)][order(chromosome_name)
]
datPosCtrlNew

