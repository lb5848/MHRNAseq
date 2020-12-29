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
#keep D0 samples
grep("D0", colnames(rowname_cts))
rowname_cts <- rowname_cts[, grep("D0", colnames(rowname_cts))]


# only integer values
rowname_cts <- round(rowname_cts)

# load and prep col.data file
coldataFile <- "coldata.txt"

colData <- read_tsv(file.path(filesPath, coldataFile))
is_tibble(colData)
names(colData)[1] <- "rowname"
col.data <- column_to_rownames(colData)

# keep D0 only
grep("D0", rownames(col.data))
col.data <- col.data[grep("D0", rownames(col.data)), ]


# 2 indicates columns, 1 indicates rows
col.data <- as.data.frame(apply(col.data, 2, as.factor))


summary(col.data)
all(colnames(rowname_cts) == rownames(col.data))

################### Running DESeq2 ##########################

dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Region)

# Pre-Filtering

dim(dds)
keep <- rowSums( counts(dds) ) >= 25
summary(keep)
dds <- dds[ keep, ]
dim(dds)

dds_tt <- dds

# varianceStabilizingTransformation
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

DESeq2::plotPCA(vsd, intgroup = "Region", ntop = 500) + stat_ellipse(type = "norm", level = 0.65)
pcaplot <- DESeq2::plotPCA(vsd, intgroup = "Region", ntop = 500) + stat_ellipse(type = "norm", level = 0.7)
class(pcaplot)
ggsave(file.path(savePath, "PCA_D0.svg"), plot = pcaplot)
ggsave(file.path(savePath, "PCA_D0.png"), plot = pcaplot)


dds$Region <- relevel(dds$Region, ref = "PBL")
dds <- DESeq(dds, test = "Wald")
resultsNames(dds)
design(dds)
results(dds)

res <- results(dds)
plotMA(res)
resLFC <- lfcShrink(dds, coef = "Region_BM_vs_PBL", type="apeglm")
plotMA(resLFC)

# ============================== tidybulk ==================================

rse <- as(dds_tt, "RangedSummarizedExperiment")

counts <- rse %>% tidybulk()
head(counts)
tail(counts)

# Plot Settings
# Use colourblind-friendly colours
friendly_cols <- dittoSeq::dittoColors()

# Set theme
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )

counts_tt <- counts %>%
  mutate(sample = str_remove(sample, "L")) %>%
  mutate(Region = str_remove(Region, "L"))

counts_filtered <- counts_tt %>% keep_abundant(factor_of_interest = Region)

counts_scaled <- counts_filtered %>% scale_abundance()


counts_scaled %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap(~ source) +
  scale_x_log10() +
  custom_theme

counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500)

counts_scal_PCA %>% pivot_sample()

# PCA plot 
counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = Region)) +
  stat_ellipse(level = 0.9) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  custom_theme

# Hierarchical clustering - heatmap 
counts_scaled %>%
  
  # extract 500 most variable genes
  keep_variable(.abundance = counts_scaled, top = 500) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = feature,
    .value = counts_scaled,
    transform = log1p
  )
