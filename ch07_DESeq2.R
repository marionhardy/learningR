
## Assumes you have run the RSubread script or
# you have a feature count file made from a .bam file
# Counts should be RAW! 

library(DESeq2)
library(purrr)
library(tidyverse)
library(stringr)
library(dplyr)

## File preparation-------------------------------------------------------------

dir <- getwd()

names <- c("C48_1","C96_1","C48_2","C96_2",
               "S48_1","S96_1","S48_2","S96_2")

files = list.files("./data/RNAseq/", pattern = ".csv")
files # check that the order of names and files is correct or you will 
# assign wrong names to samples

for(i in 1:length(files)){
  x = read_csv(paste0("./data/RNAseq/",files[i]), col_names = F)
  x = x[-1,-c(2,4:9)]
  colnames(x) = c("ensembl",paste0(names[i]))
  assign(names[i],x)
}

counts <- list(C48_1,C96_1,C48_2,C96_2,S48_1,S96_1,S48_2,S96_2) %>% 
  purrr::reduce(full_join, by = "ensembl") %>% 
  as.data.frame()

# If you get a gene version and you need to remove that info to get only
# the ENSG code

# strrep =
#   sub(pattern = "\\.(.*)","",counts$ensembl)
# 
# counts$ensembl = strrep
# counts = counts[!duplicated(counts$ensembl),]

rownames(counts) = counts$ensembl
counts = counts[,-1]
class(counts) = 'numeric'
head(counts)

# Create the coldata for the high level data structure

coldata <- data.frame(
  celltype=c(rep("ctrl",4),rep("selected",4)),
  timeline = as.factor(c(48,48,96,96,48,48,96,96)),
  replicate=rep(c(1,2),4))

rownames(coldata) <- colnames(counts)


## Expression analysis----------------------------------------------------------

# Create the DESeq object

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~celltype ) 

# Generate a linear model

dds <- DESeq(dds)


# Quality check-----------------------------------------------------------------

sizeFactors(dds) 

plotDispEsts(dds) 

resultsNames(dds) 
# dds$celltype <- relevel(dds$celltype, ref="ctrl") # donc on lui dit quelle est la ref
# dds <- DESeq(dds)
# resultsNames(dds)

# Check the distribution of residues

as_tibble(assay(dds)) %>%
  gather(sample, value = counts) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample)

# PCA preparation

rld <- rlogTransformation(dds)

plotPCA(rld,intgroup="celltype")  
plotPCA(rld,intgroup="timeline")
plotPCA(rld,intgroup="replicate") 


# Results can be extracted

res <- results(dds, name = "celltype_selected_vs_ctrl")

res_tbl <- as_tibble(res, rownames="ENSMUG")

# Save the tibbles etc

saveRDS(dds, file = "./data/210906_dds.rds")
saveRDS(res_tbl, file = "./data/210906_res_tbl.rds")
write_tsv(res_tbl, file = "./data/210906_res_tbl.tsv")

