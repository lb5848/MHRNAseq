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
workingDir <- "201228_WorkingDirectory_Final"
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

# ============================= prep for tidybulk() =======================

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

# ============================= create tidybulk() ==============================
counts <- rowname_cts %>% as_tibble(rownames = "geneID")
sampleinfo <- col.data %>% as_tibble(rownames = "sample")

counts_tt <-
  # convert to tidy format
  pivot_longer(counts, cols = !starts_with("geneID"), names_to = "sample", values_to = "counts") %>%
  
  # order the columns for tidybulk
  select(sample, geneID, counts) %>%
  
  # add the sample info
  left_join(sampleinfo) %>%
  
  # shorten sample name
  mutate(sample = str_remove(sample, "L")) %>%
  mutate(Region = str_remove(Region, "L")) %>%
  mutate(group = paste0(Region, "_", Condition)) %>%
  
  # convert to tidybulk tibble
  tidybulk(.sample = sample, .transcript = geneID, .abundance = counts)

ggplot(counts_tt, aes(x = sample, weight = counts, fill = ID)) +
  geom_bar() +
  custom_theme

# scale counts
tt_norm <- counts_tt %>% 
  keep_abundant(factor_of_interest = group, minimum_counts = 5, minimum_proportion = 0.25) %>%
  scale_abundance() 
# cannot use adjust_abundance() as only one batch is allowed

# plot scaled counts
tt_norm %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap(~ source) +
  scale_x_log10() +
  custom_theme
# boxplot
tt_norm %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = sample, y = abundance + 1, color = group)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = median(abundance + 1)), colour="red") +
  facet_wrap(~source) +
  scale_y_log10() +
  custom_theme

tt_norm_PCA <-
  tt_norm %>%
  reduce_dimensions(method = "PCA", top = 1500)

# tt_norm_PCA %>% pivot_sample()

# PCA plot 
tt_norm_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = group)) +
  stat_ellipse(level = 0.73) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  custom_theme

tt_norm_MDS <- 
  tt_norm %>%
  reduce_dimensions(method = "MDS", top = 1500)

tt_norm_MDS %>%
  pivot_sample() %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = group)) +
  stat_ellipse(level = 0.70) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  custom_theme
