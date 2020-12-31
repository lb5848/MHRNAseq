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
aMILvsaPBL <- fread(files[5]) #rev!!

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

aMILvsaPBL <- aMILvsaPBL %>%
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
  head(10)

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
  ggtitle("BM: Hypoxia vs Normoxia") + 
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("BM_HvsN_volcano.svg", plot = plot)
ggsave("BM_HvsN_volcano.png", plot = plot)

# no labels on genes!
plot <- BM_HvsN %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, "", "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  ggtitle("BM: Hypoxia vs Normoxia") + 
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("BM_HvsN_volcano_noLabels.svg", plot = plot)
ggsave("BM_HvsN_volcano_noLabels.png", plot = plot)

# =========================== PB_HvsN =======================================

topgenes <-
  PB_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(10)

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
  ggtitle("PB: Hypoxia vs Normoxia") + 
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("PB_HvsN_volcano.svg", plot = plot)
ggsave("PB_HvsN_volcano.png", plot = plot)
# no labels on genes

plot <- PB_HvsN %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, "", "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  ggtitle("PB: Hypoxia vs Normoxia") + 
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("PB_HvsN_volcano_noLabels.svg", plot = plot)
ggsave("PB_HvsN_volcano_noLabels.png", plot = plot)

# =========================== Hy_BvsP =======================================

topgenes <-
  Hy_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(10)

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
  ggtitle("Hypoxia: BM vs PB") +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("Hy_BvsP_volcano.svg", plot = plot)
ggsave("Hy_BvsP_volcano.png", plot = plot)

plot <- Hy_BvsP %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, "", "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  ggtitle("Hypoxia: BM vs PB") +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("Hy_BvsP_volcano_noLabels.svg", plot = plot)
ggsave("Hy_BvsP_volcano_noLabels.png", plot = plot)


# =========================== No_BvsP =======================================

topgenes <-
  No_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(10)

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
  ggtitle("Normoxia: BM vs PB") +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("No_BvsP_volcano.svg", plot = plot)
ggsave("No_BvsP_volcano.png", plot = plot)

plot <- No_BvsP %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, "", "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  ggtitle("Normoxia: BM vs PB") +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("No_BvsP_volcano_noLabels.svg", plot = plot)
ggsave("No_BvsP_volcano_noLabels.png", plot = plot)

# ============================= aMILs vs aPBLs ================================

topgenes <-
  aMILvsaPBL %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  arrange(padj) %>%
  head(10)

topgenes_symbols <- topgenes %>% pull(id)

plot <- aMILvsaPBL %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, as.character(id), "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  ggtitle("aMILs vs aPBLs") + 
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("aMILs_vs_aPBLs_volcano.svg", plot = plot)
ggsave("aMILs_vs_aPBLs_volcano.png", plot = plot)

# no labels on genes!
plot <- aMILvsaPBL %>%
  
  # Subset data
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  mutate(id = ifelse(id %in% topgenes_symbols, "", "")) %>%
  
  # Plot
  ggplot(aes(x = log2FoldChange, y = padj, label = id)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  ggtitle("aMILs vs aPBLs") + 
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
plot
ggsave("aMILs_vs_aPBLs_volcano_noLabels.svg", plot = plot)
ggsave("aMILs_vs_aPBLs_volcano_noLabels.png", plot = plot)

# ========================= venn ==================================
# ================== BM_HvsN PB_HvsN =============================

BMup <- BM_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(id)
dim(BMup)

PBup <- PB_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(id)

setlistup <- list( "BM Hyp vs Nor" = sample(BMup$id, length(BMup$id)), 
                   "PB Hyp vs Nor" = sample(PBup$id, length(PBup$id)))

BMdown <- BM_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(id)

PBdown <- PB_HvsN %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(id)

setlistdown <- list( "BM Hyp vs Nor" = sample(BMdown$id, length(BMdown$id)), 
                     "PB Hyp vs Nor" = sample(PBdown$id, length(PBdown$id)))

vennsetup <- overLapper(setlistup, type = "vennsets")
vennsetdown <- overLapper(setlistdown, type = "vennsets")

svg("vennBMPB_HvsN.svg")
vennPlot(list(vennsetdown, vennsetup), mymain = "DEG BM/PB Hyp vs Nor", lines = "black", lcol = "black",
         mysub = "", colmode = 2, ccol = c("blue", "red"))
dev.off()

# ================== Hyp_BvsP Norm_BvsP =============================
Hypup <- Hy_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(id)
dim(Hypup)

Norup <- No_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(id)
dim(Norup)

setlistup <- list( "Hyp BM vs PB" = sample(Hypup$id, length(Hypup$id)), 
                   "Nor BM vs PB" = sample(Norup$id, length(Norup$id)))
summary(setlistup)

Hypdown <- Hy_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(id)
dim(Hypdown)

Nordown <- No_BvsP %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(id)
dim(Nordown)

setlistdown <- list( "Hyp BM vs PB" = sample(Hypdown$id, length(Hypdown$id)), 
                   "Nor BM vs PB" = sample(Nordown$id, length(Nordown$id)))

vennsetup <- overLapper(setlistup, type = "vennsets")
vennsetdown <- overLapper(setlistdown, type = "vennsets")

svg("vennHypNor_BvsP.svg")
vennPlot(list(vennsetdown, vennsetup), mymain = "DEG Hyp/Nor BM vs PB", lines = "black", lcol = "black",
         mysub = "", colmode = 2, ccol = c("blue", "red"))
dev.off()

# =========================== enrichment/GSEA ===========================
# load msigdbr

hallmark.path <- paste(filePath, "gene_sets", "h.all.v7.1.symbols.gmt", sep = "/")
pathway.HALLMARK <- gmtPathways(hallmark.path)
HALLMARK <- read.gmt(hallmark.path)
c7.path <- paste(filePath, "gene_sets", "c7.all.v7.1.symbols.gmt", sep = "/")
pathway.C7 <- gmtPathways(c7.path)
C7 <- read.gmt(c7.path)
c5.path <- paste(filePath, "gene_sets", "c5.bp.v7.1.symbols.gmt", sep = "/")
pathway.C5 <- gmtPathways(c5.path)
C5 <- read.gmt(c5.path)
custom.path <- paste(filePath, "gene_sets", "custom.pathways.gmt", sep = "/")
pathway.CUSTOM <- gmtPathways(custom.path)
CUSTOM <- read.gmt(custom.path)

BM_HvsN_rnk <- BM_HvsN %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  dplyr::select(id, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(id) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe()

PB_HvsN_rnk <- PB_HvsN %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  dplyr::select(id, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(id) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe()

Hy_BvsP_rnk <- Hy_BvsP %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  dplyr::select(id, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(id) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe()  

No_BvsP_rnk <- No_BvsP %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  dplyr::select(id, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(id) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe() 

aMILvsaPBL_rnk <- aMILvsaPBL %>%
  mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  filter(significant) %>%
  dplyr::select(id, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(id) %>%
  summarise(stat = mean(stat)) %>%
  arrange(desc(stat)) %>%
  deframe() 

#enrich
common_Hyp <- BM_HvsN_rnk[names(BM_HvsN_rnk) %in% names(PB_HvsN_rnk)]


common_Hyp_ego <- enrichGO(names(common_Hyp), OrgDb = "org.Hs.eg.db", 
                   keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
barplot(common_Hyp_ego, showCategory = 10)

common_Hyp_gH <- GSEA(common_Hyp, exponent = 1, minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 0.25,
                      TERM2GENE = HALLMARK, by = "fgsea")
dotplot(common_Hyp_gH)
gseaplot2(common_Hyp_gH, geneSetID = "HALLMARK_HYPOXIA")

common_BM <- Hy_BvsP_rnk[names(Hy_BvsP_rnk) %in% names(No_BvsP_rnk)]

common_BM_ego <- enrichGO(names(common_BM), OrgDb = "org.Hs.eg.db", 
                          keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
barplot(common_BM_ego, showCategory = 30)

common_BM_gCustom <- GSEA(common_BM, exponent = 1, minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 0.1,
                          TERM2GENE = CUSTOM, by = "fgsea")
gseaplot2(common_BM_gCustom, geneSetID = c(3,4))
gseaplot2(common_BM_gCustom, geneSetID = "LeonePowell_UP")
gseaplot2(common_BM_gCustom, geneSetID = c("LeonePowell_DN","LeonePowell_UP"))

aMILvsaPBL_gCUSTOM <- GSEA(aMILvsaPBL_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                           TERM2GENE = CUSTOM, by = "fgsea")
aMILvsaPBL_gCUSTOM@result$ID

svg("GSEACUSTOM_aMILs_aPBLs.svg")
gseaplot2(aMILvsaPBL_gCUSTOM, geneSetID = c("IL21LDHivsIL21_UP", "IL2LDHivsIL2_UP"), 
          title = "aMILs vs aPBLs", ES_geom = "line")
dev.off()

aMILvsaPBL_gH <- GSEA(aMILvsaPBL_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                           TERM2GENE = HALLMARK, by = "fgsea")
aMILvsaPBL_gH@result$ID
svg("GSEAHALLMARK_aMILs_aPBLs.svg")
gseaplot2(aMILvsaPBL_gH, geneSetID = c("HALLMARK_HYPOXIA" , "HALLMARK_MTORC1_SIGNALING", 
                                       "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_P53_PATHWAY"  ), 
          title = "aMILs vs aPBLs", ES_geom = "line")
dev.off()

aMILvsaPBL_gC5 <- GSEA(aMILvsaPBL_rnk, exponent = 1, minGSSize = 15, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                           TERM2GENE = C5, by = "fgsea")
aMILvsaPBL_gC5@result$ID
aMILvsaPBL_gC5@result$ID[grep("RESPONSE_TO", aMILvsaPBL_gC5@result$ID)]

svg("GSEAC5_aMILs_aPBLs.svg")
gseaplot2(aMILvsaPBL_gC5, geneSetID = c("GO_RESPONSE_TO_CYTOKINE", "GO_RESPONSE_TO_INTERFERON_GAMMA" ), 
          title = "aMILs vs aPBLs", ES_geom = "line")
dev.off()

aMILvsaPBL_gC7 <- GSEA(aMILvsaPBL_rnk, exponent = 1, minGSSize = 15, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                       TERM2GENE = C7, by = "fgsea")
head(aMILvsaPBL_gC7@result$ID, n = 10)
aMILvsaPBL_gC7@result$ID[grep("EXHA", aMILvsaPBL_gC7@result$ID)]

svg("GSEAC7_aMILs_aPBLs_EXHAUSTION.svg")
gseaplot2(aMILvsaPBL_gC7, geneSetID = c(
                                        "GSE41867_DAY15_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_UP",
                                        "GSE41867_DAY15_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_DN",
                                        "GSE41867_DAY6_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_DN",
                                        "GSE41867_DAY6_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_UP"),
          title = "aMILs vs aPBLs", ES_geom = "line")
dev.off()

svg("GSEAC7_aMILs_aPBLs_STEM.svg")
gseaplot2(aMILvsaPBL_gC7, geneSetID = c("GSE41867_MEMORY_VS_EXHAUSTED_CD8_TCELL_DAY30_LCMV_UP", 
                                        "GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_DN",
                                        "GSE23321_CENTRAL_MEMORY_VS_NAIVE_CD8_TCELL_UP"), 
          title = "aMILs vs aPBLs", ES_geom = "line")
dev.off()

Hyp_BvP_gCUSTOM <- GSEA(Hy_BvsP_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
              TERM2GENE = CUSTOM, by = "fgsea")
Hyp_BvP_gCUSTOM@result$ID
gseaplot2(Hyp_BvP_gCUSTOM, geneSetID =c("LeonePowell_UP", "LeonePowell_DN"), ES_geom = "dot")

svg("GSEACUSTOM_BMH_PBH.svg")
gseaplot2(Hyp_BvP_gCUSTOM, geneSetID =c("LeonePowell_DN", "LeonePowell_UP", "IL2LDHivsIL2_UP"),
          title = "GSEA: BM_Hyp vs PB_Hyp - CUSTOM gmt", pvalue_table = FALSE, ES_geom = "dot")
dev.off()

Hy_BvP_gC5 <- GSEA(Hy_BvsP_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                   TERM2GENE = C5, by = "fgsea")
Hy_BvP_gC5@result$ID[grep("APOPT", Hy_BvP_gC5@result$ID)]
Hy_BvP_gC5@result$ID

svg("GOBP_HypBMvsPB_Tcell.svg")
gseaplot2(Hy_BvP_gC5, 
          geneSetID =c("GO_CYTOKINE_PRODUCTION", 
                       "GO_T_CELL_DIFFERENTIATION",
                       "GO_CHEMOKINE_PRODUCTION"), 
          ES_geom = "line", title = "GSEA: Hypoxia - MILs vs PBLs - GOBP")
dev.off()
svg("GOBP_HypBMvsPB_APOPTOSIS.svg")
gseaplot2(Hy_BvP_gC5, 
          geneSetID =c("GO_APOPTOTIC_SIGNALING_PATHWAY", 
                       "GO_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
                       "GO_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY" ), 
          ES_geom = "line", title = "GSEA: Hypoxia - MILs vs PBLs - GOBP")
dev.off()

Hy_BvP_gC7 <- GSEA(Hy_BvsP_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                   TERM2GENE = C7, by = "fgsea")
Hy_BvP_gC7@result$ID[grep("EXHAUSTED", Hy_BvP_gC7@result$ID)]

svg("GSEAC7_Hyp_BvsP.svg")
gseaplot2(Hy_BvP_gC7, 
          geneSetID =c("GSE23321_CD8_STEM_CELL_MEMORY_VS_EFFECTOR_MEMORY_CD8_TCELL_DN",
                       "GSE41867_DAY6_EFFECTOR_VS_DAY30_EXHAUSTED_CD8_TCELL_LCMV_CLONE13_DN",
                       "GSE41867_MEMORY_VS_EXHAUSTED_CD8_TCELL_DAY30_LCMV_UP"), 
          ES_geom = "line", title = "GSEA: Hypoxia - MILs vs PBLs - C7")
dev.off()

Hy_BvP_gH <- GSEA(Hy_BvsP_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                   TERM2GENE = HALLMARK, by = "fgsea")
Hy_BvP_gH@result$ID

BM_HvN_gCUSTOM <- GSEA(BM_HvsN_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                          TERM2GENE = CUSTOM, by = "fgsea")
BM_HvN_gCUSTOM@result$ID

svg("GSEACUSTOM_BM_HvsN.svg")
gseaplot2(Hyp_BvP_gCUSTOM, geneSetID =c("LeonePowell_UP", "IL2vsIL21_UP"), ES_geom = "dot")
dev.off()

BM_HvN_gC5 <- GSEA(BM_HvsN_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                       TERM2GENE = C5, by = "fgsea")
BM_HvN_gC5@result$ID

BM_HvN_gC7 <- GSEA(BM_HvsN_rnk, exponent = 1, minGSSize = 5, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                   TERM2GENE = C7, by = "fgsea")
BM_HvN_gC7@result$ID