rm(list = ls())

# load libraries

library(DESeq2)

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
# =======================Prep data for tidybulk================================
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

# prep DESeq2 object 

dds <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                              design = ~ ID + Condition)
# setup multifactorial design

# create "group" - ?levels "BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"

dds$group <- factor(paste0(dds$Region, "_", dds$Condition),
                    levels = c("PBL_Norm", "BM_Norm", "PBL_Hyp", "BM_Hyp"))
design(dds) <- formula(~ group + ID)


# coerce DESeqDataSet to RangedSummarizedExperiment
rse <- as(dds, "RangedSummarizedExperiment")

counts <- rse %>% tidybulk()
head(counts)
tail(counts)

# remove "L" from "PBL"
counts_format <- counts %>% 
  mutate(Region = str_remove(Region, "L")) %>%
  mutate(group = str_remove(group, "L"))
levels(factor(counts_format$group))
levels(factor(counts_format$Region))

counts_scaled <- counts_format %>%
  identify_abundant(factor_of_interest = group, minimum_counts = 25, minimum_proportion = 0.25) %>%
  scale_abundance(method = "TMM")

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
  ggplot(aes(x = PC1, y = PC2, colour = group, shape = Region)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample), show.legend = FALSE) +
  theme_bw()

# Reduce data dimensionality with arbitrary number of dimensions
tt_mds <- counts_scaled %>% reduce_dimensions(method = "MDS", .dims = 6, top = 1500)

tt_mds %>%
  pivot_sample() %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = group, shape = Region)) +
  geom_point(size = 4) +
  stat_ellipse(level = 0.65, type = "norm") +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  theme_bw()

counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract 500 most variable genes
  keep_variable( .abundance = counts_scaled, top = 40) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = feature,
    .value = counts_scaled,
    annotation = c(group),
    transform = log1p
  )

counts_de <- counts_scaled %>%
  test_differential_abundance(
    .formula = ~ 0 + group + ID,
    .contrasts = c("groupBM_Hyp - groupPB_Hyp"),
    omit_contrast_in_colnames = TRUE
  )

counts_de %>% pivot_transcript(.transcript = feature)

topgenes <-
  counts_de %>%
  pivot_transcript() %>%
  arrange(PValue) %>%
  head(20)

topgenes_symbols <- topgenes %>% pull(feature)

counts_de %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
