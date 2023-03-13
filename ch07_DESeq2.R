
## Assumes you have run the RSubread script or
# you have a feature count file made from a .bam file
# Counts should be RAW! 

library(DESeq2)
library(purrr)
library("tidyverse")
library(stringr)
library(dplyr)

## File preparation-------------------------------------------------------------

dir <- getwd()

C48_1 <- read_csv(paste0(dir,"/data/counts/Ctrl_1_48h_counts.csv"))
C48_2 <- read_csv(paste0(dir,"/data/counts/Ctrl_2_48h_counts.csv"))
C96_1 <- read_csv(paste0(dir,"/data/counts/Ctrl_1_96h_counts.csv"))
C96_2 <- read_csv(paste0(dir,"/data/counts/Ctrl_2_96h_counts.csv"))
S48_1 <- read_csv(paste0(dir,"/data/counts/Sel_1_48h_counts.csv"))
S48_2 <- read_csv(paste0(dir,"/data/counts/Sel_2_48h_counts.csv"))
S96_1 <- read_csv(paste0(dir,"/data/counts/Sel_1_96h_counts.csv"))
S96_2 <- read_csv(paste0(dir,"/data/counts/Sel_2_96h_counts.csv"))

datanames <- c("C48_1","C48_2","C96_1","C96_2",
               "S48_1","S48_2","S96_1","S96_2")

list <- llist(D0_1,D0_2,C48_1,C48_2,C96_1,C96_2,S48_1,S48_2,S96_1,S96_2)

# Make count matrix (keep gene_id as rownames and rawcounts)

fun <- function(i){
  i <- i %>% as.data.frame()
  rownames(i) <- i$ensembl_gene_id
  i <- i %>% dplyr::select(RawCounts, ensembl_gene_id)
}

list <- lapply(list,fun)
list2env(list, envir=.GlobalEnv)

head(D0_1)

# Change rawcounts to sample name

colnames(C48_1)[1] <- datanames[3]
colnames(C48_2)[1] <- datanames[4]
colnames(C96_1)[1] <- datanames[5]
colnames(C96_2)[1] <- datanames[6]
colnames(S48_1)[1] <- datanames[7]
colnames(S48_2)[1] <- datanames[8]
colnames(S96_1)[1] <- datanames[9]
colnames(S96_2)[1] <- datanames[10]

counts <- list(C48_1,C48_2,C96_1,C96_2,S48_1,S48_2,S96_1,S96_2) %>% 
  purrr::reduce(dplyr::full_join,)

rownames(counts) <- counts$X
counts <- counts %>% dplyr::select(-X)
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
  ggplot(aes(x = counts, fill = sample)) +
  geom_histogram(bins = 20) 
facet_wrap(~ sample)

as_tibble(assay(dds)) %>%
  gather(sample, value = counts) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample)

# PCA preparation

rld <- rlogTransformation(dds)

plotPCA(rld,intgroup="celltype")  
plotPCA(rld,intgroup="timepoint")
plotPCA(rld,intgroup="replicate") 


# Results can be extracted

res <- results(dds, name = "celltype_selected_vs_ctrl")

res_tbl <- as_tibble(res, rownames="ENSMUG")

# Save the tibbles etc

saveRDS(dds, file = "./data/210906_dds.rds")
saveRDS(res_tbl, file = "./data/210906_res_tbl.rds")
write_tsv(res_tbl, file = "./data/210906_res_tbl.tsv")

