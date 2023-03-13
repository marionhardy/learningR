##############
## Q and A ##
#############

library(tidyverse)
library(ggplot2)

# Plotting a volcano plot-------------------------------------------------------

## Load the data we will use for this course
## RNAseq data for volcano plot + highlight per filtering

rna <- read_csv("./data/res_tbl.csv")
str(rna)
head(rna)

# To plot RNAseq data, a violinplot is very common

rna %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 0.5) 

rna %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 0.5)+
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  labs(title = "Selected cells T96 vs T48")+
  theme_bw()

rna %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
              color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  labs(title = "Selected cells T96 vs T48")+
  theme_bw()

rna %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|padj<0.05&log2FoldChange<=-1,
                            as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  geom_text_repel(max.overlaps = 10)+
  labs(title = "Selected cells T96 vs T48")+
  theme_bw()


## For your information, base R also has plot functions which function
# very well and can be used when doing quality control where there is no
# need for fancy colors and titles


hist(rna$padj)

plot(rna$log2FoldChange, -log10(rna$padj))
abline(a = -log10(0.05), b= 0, v = c(1,-1))

boxplot(airquality$Temp ~ airquality$Measurer)

dev.null()


# High-level data structures----------------------------------------------------

# Besides lists, dataframesn vectors... You will encounter other types of
# objects while using R, especially in using R for omics analyses

# Some of these objects include:

# Rectangular feature x sample data – 
# SummarizedExperiment::SummarizedExperiment() (RNAseq count matrix, microarray, …)
# 
# Genomic coordinates – GenomicRanges::GRanges() (1-based, closed interval)
# DNA / RNA / AA sequences – Biostrings::*StringSet()
# Multi-omics data – MultiAssayExperiment::MultiAssayExperiment()
# Single cell data – SingleCellExperiment::SingleCellExperiment()
# Quantitative proteomics data – QFeatures::QFeatures()

# Tidyverse cannot be used on those objects, subsetting is made by using the
# brackets

BiocManager::install("SummarizedExperiment")
BiocManager::install("airway")

airway_df <- read_csv("./data/ch04_airway.csv")
























