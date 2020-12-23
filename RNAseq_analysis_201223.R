rm(list = ls())

# load libraries

library(DESeq2)
library(hciR)

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




# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create working directory
workingDir <- "201220_WorkingDirectory_tidytranscriptomics"
dirPath <- file.path(PrimaryDirectory, workingDir)
dir.create(dirPath)
setwd(dirPath)

savePath <- file.path(dirPath, "results")
saveDir <- file.path(dirPath, "results")
dir.create(savePath)

# load files
filesPath <- file.path(dirPath, "files")
countFile <- "Counts"

# load and prep cts file
cts <- read_tsv(file.path(filesPath, countFile))
is_tibble(cts)
rowname_cts <- cts
names(rowname_cts)[1] <- "rowname"
rowname_cts <- column_to_rownames(rowname_cts)
colnames(rowname_cts)
# remove D0 samples
grep("D0", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("D0", colnames(rowname_cts))]

# remove DA samples
grep("DA", colnames(rowname_cts))
rowname_cts <- rowname_cts[, -grep("DA", colnames(rowname_cts))]
colnames(rowname_cts)


# only integer values
rowname_cts <- round(rowname_cts)

# load and prep col.data file
coldataFile <- "coldata.txt"

colData <- read_tsv(file.path(filesPath, coldataFile))
is_tibble(colData)
names(colData)[1] <- "rowname"
col.data <- column_to_rownames(colData)

# remove D0 and DA samples
grep("D0", rownames(col.data))
col.data <- col.data[-grep("D0", rownames(col.data)), ]

grep("DA", rownames(col.data))
col.data <- col.data[-grep("DA", rownames(col.data)), ]
rownames(col.data)


# 2 indicates columns, 1 indicates rows
col.data <- as.data.frame(apply(col.data, 2, as.factor))


summary(col.data)
all(colnames(rowname_cts) == rownames(col.data))

################### Running DESeq2 ##########################

dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Condition)
# setup multifactorial design

# create "group" - ?levels "BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"
dds$group <- factor(paste0(dds$Region, "_", dds$Condition),
                    levels = c("BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"))
design(dds) <- formula(~ ID + group)

# Pre-Filtering
dim(dds)
keep <- rowSums( counts(dds) ) >= 25
summary(keep)
dds <- dds[ keep, ]
dim(dds)

# varianceStabilizingTransformation
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

DESeq2::plotPCA(vsd, intgroup = "group", ntop = 500)
ggsave(file.path(savePath, "PCA_4pts_activated.svg"), plot = last_plot())
ggsave(file.path(savePath, "PCA_4pts_activated.png"), plot = last_plot())

# ============================== tidybulk ==================================
# coerce DESeqDataSet to RangedSummarizedExperiment
rse <- as(dds, "RangedSummarizedExperiment")

counts <- rse %>% tidybulk()
head(counts)
tail(counts)

# # remove "L" from "PBL"
# counts_format <- counts %>% 
#   mutate(Region = str_remove(Region, "L")) %>%
#   mutate(group = str_remove(group, "L"))
# levels(factor(counts_format$group))
# levels(factor(counts_format$Region))

counts_scaled <- counts_format %>%
  identify_abundant(factor_of_interest = group, minimum_counts = 10, minimum_proportion = 0.25) %>%
  keep_abundant() %>% keep_variable(top = 500) %>%
  scale_abundance()


counts_scaled %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap(~source) +
  scale_x_log10() +
  theme_bw()

counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500)

counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  theme_bw()

counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract 500 most variable genes
  keep_variable( .abundance = counts_scaled, top = 500) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = feature,
    .value = counts_scaled,
    annotation = c(group),
    transform = log1p
  )

############################ back to DESeq2 ################################

dds <- DESeq(dds, test = "Wald")
resultsNames(dds)
design(dds)

res <- results_all(dds, vs = "all", trt = "group")

p.adj.cutoff <- 0.01
log2FC.cutoff <- 2
res_sig <- res %>% lapply(function(x) {
  x <- x %>% as_tibble() %>% filter(padj < p.adj.cutoff) %>% 
    filter(abs(log2FoldChange) >= log2FC.cutoff)})
res_sig

# save csv files - significant genes only - all contrasts
csvPath <- file.path(savePath, "csv_files")
dir.create(csvPath)
lapply(1:length(res_sig), function(i){
  res_sig[[i]] %>% as_tibble() %>% arrange(padj) %>% 
    fwrite(file.path(csvPath, paste0(names(res_sig[i]), ".csv")))
})
write_deseq(res_sig, dds, vsd, file = file.path(csvPath, "DESeq2_results.xlsx"))
