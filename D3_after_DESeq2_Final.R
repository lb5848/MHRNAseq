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

setwd(savePath)

# load file

csvPath <- file.path(savePath, "csv_files")
setwd(csvPath)

files <- list.files(path = ".", pattern = ".csv$", full.names = FALSE)
files <- files[grep("all_", files)]
Hy_BvsP <- fread(files[1])
BM_HvsN <- fread(files[2]) # rev!!
No_BvsP <- fread(files[4])
PB_HvsN <- fread(files[6]) # rev!!

plotPath <- file.path(csvPath, "plots")
dir.create(plotPath)
setwd(plotPath)

# rev!!
BM_HvsN <- BM_HvsN %>% 
  mutate(log2FoldChange = (-1) * log2FoldChange) %>%
  mutate(stat = (-1) * stat)

PB_HvsN <- PB_HvsN %>%
  mutate(log2FoldChange = (-1) * log2FoldChange) %>%
  mutate(stat = (-1) * stat)

# ============================ Plot Settings ==============================
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


# make volcano plots

# ========================== volcano ==================================
 # ============================ BM_HvsN ==================================
topgenes <-
  BM_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(20)

topgenes_symbols <- topgenes %>% pull(id)

plot <- BM_HvsN %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, as.character(id), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("BM_HvsN_volcano.svg", plot = plot)
ggsave("BM_HvsN_volcano.png", plot = plot)
# =========================== PB_HvsN =======================================

topgenes <-
  PB_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(20)

topgenes_symbols <- topgenes %>% pull(id)

plot <- PB_HvsN %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, as.character(id), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("PB_HvsN_volcano.svg", plot = plot)
ggsave("PB_HvsN_volcano.png", plot = plot)

# =========================== Hy_BvsP =======================================

topgenes <-
  Hy_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(20)

topgenes_symbols <- topgenes %>% pull(id)

plot <- Hy_BvsP %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, as.character(id), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("Hy_BvsP_volcano.svg", plot = plot)
ggsave("Hy_BvsP_volcano.png", plot = plot)


# =========================== No_BvsP =======================================

topgenes <-
  No_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(20)

topgenes_symbols <- topgenes %>% pull(id)

plot <- No_BvsP %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, as.character(id), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("No_BvsP_volcano.svg", plot = plot)
ggsave("No_BvsP_volcano.png", plot = plot)

# ========================= venn ==================================
# ================== BM_HvsN PB_HvsN =============================
BMgenes <- BM_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  select(id)
dim(BMgenes)


PBgenes <- PB_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  select(id)
dim(PBgenes)

comb <- c(BMgenes$id, PBgenes$id)
res_BM_HvsN <- comb %in% BMgenes$id
res_PB_HvsN <- comb %in% PBgenes$id

res <- vennCounts(cbind(res_BM_HvsN, res_PB_HvsN))

png("vennBMPB_HvsN.png")
vennDiagram(res, cex = 1, names = c("BM Hyp vs Norm","PB Hyp vs Norm"), circle.col = c("red", "blue"))
dev.off()

# ================== Hyp_BvsP Norm_BvsP =============================
Hypgenes <- Hy_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  select(id)
dim(Hypgenes)


Norgenes <- No_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  select(id)
dim(Norgenes)

comb <- rbind(Hypgenes, Norgenes)
res_Hyp <- comb$id %in% Hypgenes$id
res_Nor <- comb$id %in% Norgenes$id


res <- vennCounts(cbind(res_Hyp, res_Nor))
png("vennHyp_BvsP.png")
vennDiagram(res, cex = 1, names = c("Hyp BM vs PB","Nor BM vs PB"), circle.col = c("red", "blue"))
dev.off()



# =========================== enrichment/GSEA ===========================