rm(list = ls())

install_packages <- FALSE

if(install_packages){
  BiocManager::install("DOSE")
  BiocManager::install("clusterProfiler")
  BiocManager::install("enrichplot")
  BiocManager::install("ggupset")
  BiocManager::install("gmt")
  BiocManager::install("hypeR")
  BiocManager::install("gprofiler2")
}

library(DOSE)
library(clusterProfiler)
library(tidyverse)
library(readr)
library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(hypeR)
library(gprofiler2)
library(magrittr)
library(data.table)
library(biomaRt)
library(readxl)
library(enrichplot)
library(org.Hs.eg.db)
library(ggupset)
library(gmt)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create working directory
workingDir <- "201220_WorkingDirectory_tidytranscriptomics"
dirPath <- file.path(PrimaryDirectory, workingDir)
setwd(dirPath)

resPath <- file.path(dirPath, "results")
resDir <- file.path(dirPath, "results")
filePath <- file.path(resPath, "csv_files")

setwd(filePath)

gseaPath <- file.path(filePath, "enrichment_results")
dir.create(gseaPath)

files <- list.files(filePath, pattern = ".csv$", full = FALSE)

BMHvsBMN <- fread(files[1])
rnk <- BMHvsBMN %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(id, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(id) %>%
  summarise(stat = mean(stat))
ranks <- rnk %>% arrange(desc(stat))
# ranks$stat <- rnk$stat * -1
ranks <- deframe(ranks)

kk <- gseKEGG(ranks, organism = "hsa")
# convert mouseGenes to human SYMBOL
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
c7.path <- paste(filePath, "gene_sets", "c7.all.v7.1.symbols.gmt", sep = "/")
pathway.C7 <- gmtPathways(c7.path)
head(pathway.C7)
c5.path <- paste(filePath, "gene_sets", "c5.bp.v7.1.symbols.gmt", sep = "/")
pathway.C5 <- gmtPathways(c5.path)
  

fgseaRes <- fgsea(pathway.C5, ranks, maxSize = 500, minSize = 20, eps = 0)
# plotEnrichment(pathway.HALLMARK[["HALLMARK_HYPOXIA"]], ranks) + labs(title = "HALLMARK_HYPOXIA")

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(padj < 0.01) %>%
  filter(abs(NES) >= 2) %>%
  filter(!grepl("MYELOID", pathway))


fgseaResTidy

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
