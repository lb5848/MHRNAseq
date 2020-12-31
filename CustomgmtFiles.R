rm(list = ls())

# load libraries

library(DESeq2)
library(hciR)
library(ggplotify)
library(data.table)

# tidyverse core packages
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)

# tidyverse-friendly packages
library(tidyHeatmap)
library(tidybulk)
library(ggrepel)
library(plotly)
library(GGally)

# colorblind-friendly packages
library(dittoSeq)

# venn diagram packages
library(VennDiagram)
library(systemPipeR)
library(ggVennDiagram)
library(limma)

# gene enrichment analysis packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(hciRdata)
library(enrichplot)
library(DOSE)
library(fgsea)
library(biomaRt)

# manage xls files
library(readxl)



# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create working directory
workingDir <- "201228_WorkingDirectory_Final"
dirPath <- file.path(PrimaryDirectory, workingDir)
dir.create(dirPath)
setwd(dirPath)

savePath <- file.path(dirPath, "results")
saveDir <- file.path(dirPath, "results")
dir.create(savePath)

resPath <- file.path(dirPath, "results")
resDir <- file.path(dirPath, "results")
filePath <- file.path(resPath, "csv_files")


setwd(filePath)

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                   values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = TRUE)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


hallmark.path <- paste(filePath, "gene_sets", "h.all.v7.1.symbols.gmt", sep = "/")
pathway.HALLMARK <- gmtPathways(hallmark.path)
HALLMARK <- read.gmt(hallmark.path)

c7.path <- paste(filePath, "gene_sets", "c7.all.v7.1.symbols.gmt", sep = "/")
pathway.C7 <- gmtPathways(c7.path)
C7 <- read.gmt(c7.path)
head(pathway.C7)
c5.path <- paste(filePath, "gene_sets", "c5.bp.v7.1.symbols.gmt", sep = "/")
pathway.C5 <- gmtPathways(c5.path)
C5 <- read.gmt(c5.path)

powellXls <- read_excel(paste(filePath, "gene_sets", "aav2588_Leone_DataS1.xlsx", sep = "/"), 
                        sheet = 2, col_names = TRUE)
colnames(powellXls)[1] <- "Gene"

PowellUP <- powellXls %>% filter(`Log2(fold_change)` > 0) %>% dplyr::select(Gene)
PowellDN <- powellXls %>% filter(`Log2(fold_change)` < 0) %>% dplyr::select(Gene)

hPowellUP <- PowellUP$Gene %>% convertMouseGeneList()
hPowellDN <- PowellDN$Gene %>% convertMouseGeneList()

xls.pnas <- "pnas.1920413117.sd01.xlsx"
xls.pnas <- paste(filePath, "gene_sets", xls.pnas, sep = "/")

IL2vsNC_UP <- read_excel(xls.pnas, sheet = 2, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL2vsNC_UP <- IL2vsNC_UP$Symbol
IL2vsNC_UP <- convertMouseGeneList(IL2vsNC_UP)

IL21vsNC_UP <- read_excel(xls.pnas, sheet = 3, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL21vsNC_UP <- IL21vsNC_UP$Symbol
IL21vsNC_UP <- convertMouseGeneList(IL21vsNC_UP)

IL2LDHivsIL2_UP <- read_excel(xls.pnas, sheet = 4, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL2LDHivsIL2_UP <- IL2LDHivsIL2_UP$Symbol
IL2LDHivsIL2_UP <- convertMouseGeneList(IL2LDHivsIL2_UP)

IL21LDHivsIL21_UP <- read_excel(xls.pnas, sheet = 5, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL21LDHivsIL21_UP <- IL21LDHivsIL21_UP$Symbol
IL21LDHivsIL21_UP <- convertMouseGeneList(IL21LDHivsIL21_UP)

IL2vsIL21_UP <- read_excel(xls.pnas, sheet = 6, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL2vsIL21_UP <- IL2vsIL21_UP$Symbol
IL2vsIL21_UP <- convertMouseGeneList(IL2vsIL21_UP)

IL2LDHivsIL21LDHi_UP <- read_excel(xls.pnas, sheet = 7, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL2LDHivsIL21LDHi_UP <- IL2LDHivsIL21LDHi_UP$Symbol
IL2LDHivsIL21LDHi_UP <- convertMouseGeneList(IL2LDHivsIL21LDHi_UP)

GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP <- pathway.C7[["GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP"]]

varNames <- ls()[grep("_UP", ls())]
varNames <- c(varNames, "LeonePowell_UP", "LeonePowell_DN")

length(varNames)
customPathways <- pathway.HALLMARK[1:9]
names(customPathways) <- varNames
names(customPathways)
customPathways$GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP <- GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP
customPathways$IL21LDHivsIL21_UP <- IL21LDHivsIL21_UP
customPathways$IL21vsNC_UP <- IL21vsNC_UP
customPathways$IL2LDHivsIL2_UP <- IL2LDHivsIL2_UP
customPathways$IL2LDHivsIL21LDHi_UP <- IL2LDHivsIL21LDHi_UP
customPathways$IL2vsIL21_UP <- IL2vsIL21_UP
customPathways$IL2vsNC_UP <- IL2vsNC_UP
customPathways$LeonePowell_UP <- hPowellUP
customPathways$LeonePowell_DN <- hPowellDN

writeGmtPathways(customPathways, paste(filePath, "gene_sets", "custom.pathways.gmt", sep = "/"))
customPathways <- read.gmt(paste(filePath, "gene_sets", "custom.pathways.gmt", sep = "/"))
