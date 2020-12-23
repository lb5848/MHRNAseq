rm(list = ls())

# load libraries

library(DESeq2)
library(hciR)
library(ggplotify)

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
pcaplot <- DESeq2::plotPCA(vsd, intgroup = "group", ntop = 500)
class(pcaplot)
ggsave(file.path(savePath, "PCA_4pts_activated.svg"), plot = pcaplot)
ggsave(file.path(savePath, "PCA_4pts_activated.png"), plot = pcaplot)

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

ddsLRT <- dds
ddsLRT <- DESeq(ddsLRT, test = "LRT", reduced = ~ ID)

resultsNames(ddsLRT)
resLRT <- results(ddsLRT, name = resultsNames(ddsLRT)[6])
p.adj.cutoff <- 0.01
log2FC.cutoff <- 2

sig_resLRT <- resLRT %>% data.frame() %>% 
  rownames_to_column(var = "id") %>%
  as_tibble() %>%
  filter(padj < p.adj.cutoff) %>%
  filter(abs(log2FoldChange) >= log2FC.cutoff)

sig_resLRThciR <- sig_resLRT %>% data.frame() %>% column_to_rownames( var = "id" ) %>%
  rownames_to_column(var = "id") %>%
  as_tibble()

x <- top_counts(sig_resLRThciR, vsd, top = 400, filter = TRUE, sort_fc = TRUE)

plot_genes(x, intgroup = "group", scale = "diff", show_rownames = FALSE, annotation_names_col = FALSE, 
           show_colnames = TRUE)
plot <- plot_genes(x, intgroup = "group", scale = "diff", show_rownames = FALSE, 
                   annotation_names_col = FALSE, show_colnames = FALSE, output = "pheatmap")
plot <- as.ggplot(plot, scale = 1, hjust = 0, vjust = 0)
ggsave(file.path(savePath, "heatmap400genesBMHvsN.svg"), plot = plot)
ggsave(file.path(savePath, "heatmap400genesBMHvsN.png"), plot = plot)
# ============================== tidybulk ==================================
# coerce DESeqDataSet to RangedSummarizedExperiment
dds_tt <- DESeqDataSetFromMatrix(countData = rowname_cts, colData = col.data, 
                                 design = ~ ID + Condition)
# setup multifactorial design

# create "group" - ?levels "BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"
dds_tt$group <- factor(paste0(dds_tt$Region, "_", dds_tt$Condition),
                    levels = c("BM_Norm", "PBL_Norm", "BM_Hyp", "PBL_Hyp"))
design(dds_tt) <- formula(~ ID + group)
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
  mutate(Region = str_remove(Region, "L")) %>%
  mutate(group = str_remove(group, "L"))
  
counts_filtered <- counts_tt %>% keep_abundant(factor_of_interest = group)

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
  reduce_dimensions(method = "PCA", top = 1500)

counts_scal_PCA %>% pivot_sample()

# PCA plot 
counts_scal_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  custom_theme


