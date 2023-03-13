
library(DESeq2)
library(GEOquery)
library(tidyverse)

## Getting the files from NIH.gov ----------------------------------------------

gse <-  getGEO("GSE189685", GSEMatrix = F)
head(Meta(gse))
names(GSMList(gse))
GSMList(gse)[[1]]

gse189685 <- getGEO('GSE189685',GSEMatrix=TRUE)
show(gse189685)
show(pData(phenoData(gse189685[[1]]))[1:9,c(1,6,8)])

gsm <- gse189685[[1]]
exprs(gsm) # We are missing the count data
pData(gsm) # It's available as a .csv supplementary file for each GSMxxxxxxx

# Extract everything in the same folder-----------------------------------------

supp_files <- getGEOSuppFiles(GEO = "GSE189685", 
                              baseDir = paste0(getwd(),"/data")) 
# or download them manually

gunzip("./data/GSE189685/GSE189685_Raw_gene_counts.tsv.gz", overwrite = F, remove = FALSE)

## Getting the counts matrix----------------------------------------------------

counts <- read.table("./data/GSE189685/GSE189685_Raw_gene_counts.tsv", 
                     row.names = 1)

###############
## EXERCISES ##
###############

## 2----------------------------------------------------------------------------

head(counts)
colnames(counts)
rownames(counts)

coln <- counts[1,]
colnames(counts) <- coln

# Or

coln <- c("PBS_1","PBS_2","PBS_3",
          "231_1","231_2","231_3",
          "NME1_1","NME1_2","NME1_3")
colnames(counts) <- coln

head(counts)

counts <- counts[-1,]

head(counts)

str(counts)

counts <-  type.convert(counts, as.is = TRUE)

str(counts)

## 3----------------------------------------------------------------------------

coldata <- data.frame(
  celltype=rep("MCF7",9),
  treatment = c(rep("PBS",3),rep("231",3), rep("NME1",3)),
  timepoint = as.factor(rep(24,9)),
  replicate=as.factor(rep(c(1:3),3)))

rownames(coldata) <- colnames(counts)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~treatment) 


dds <- DESeq(dds)


# 4 & 5-------------------------------------------------------------------------

sizeFactors(dds) 

plotDispEsts(dds) 

resultsNames(dds) 

dds$treatment <- relevel(dds$treatment, ref="PBS")
dds <- DESeq(dds)
resultsNames(dds)

sizeFactors(dds) 

plotDispEsts(dds) 

resultsNames(dds) 

# 6-----------------------------------------------------------------------------

rld <- rlogTransformation(dds)

plotPCA(rld,intgroup="celltype")  
plotPCA(rld,intgroup="timepoint")

p1 <- plotPCA(rld,intgroup="replicate") 
p2 <- plotPCA(rld,intgroup="treatment") 

ggsave("PCA_replicate.png", p1, "./figures", width = 5, height = 5,
       dpi = 300, device = png)
ggsave("PCA_treatment.png", p2, "./figures", width = 5, height = 5,
       dpi = 300, device = png)


# 7-----------------------------------------------------------------------------

res <- results(dds, name = "treatment_NME1_vs_PBS")

res_tbl <- as_tibble(res, rownames="Ensembl")

# Save the tibbles etc

saveRDS(dds, file = "./data/210906_dds.rds")
saveRDS(res_tbl, file = "./data/210906_res_tbl.rds")
write_tsv(res_tbl, file = "./data/210906_res_tbl.tsv")


as_tibble(counts(dds["PRKAA2"], normalize = TRUE),
          rownames = 'Ensembl') %>%
  gather(sample, counts, -Ensembl) %>%
  left_join(as_tibble(coldata, rownames = "sample")) %>%
  ggplot(aes(x = sample, y = counts, fill = treatment)) +
  geom_bar(stat = 'identity', color = "grey30", width = .75) +
  theme(axis.text.x = element_text(size = 11, angle = 1),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11))+
  ggtitle("PRKAA2 mRNA counts")

as_tibble(counts(dds["PRKAA1"], normalize = TRUE),
          rownames = 'Ensembl') %>%
  gather(sample, counts, -Ensembl) %>%
  left_join(as_tibble(coldata, rownames = "sample")) %>%
  ggplot(aes(x = sample, y = counts, fill = treatment)) +
  geom_bar(stat = 'identity', color = "gray30", width = .75) +
  theme(axis.text.x = element_text(size = 11, angle = 1),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11))+
  ggtitle("PRKAA1 mRNA counts")















































