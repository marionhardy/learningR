###############################
## ch08 RNAseq data analysis ##
###############################

BiocManager::install("DESeq2")
BiocManager::install("purrr")
BiocManager::install("stringr")
BiocManager::install("Hmisc")

library(DESeq2)
library(purrr)
library(tidyverse)
library(stringr)
library(dplyr)
library(Hmisc)

## File preparation-------------------------------------------------------------

dir <- getwd()

C96_1 <- read_csv(paste0(dir,"/data/counts/Ctrl_1_96h_counts.csv"))
C96_2 <- read_csv(paste0(dir,"/data/counts/Ctrl_2_96h_counts.csv"))
S96_1 <- read_csv(paste0(dir,"/data/counts/Sel_1_96h_counts.csv"))
S96_2 <- read_csv(paste0(dir,"/data/counts/Sel_2_96h_counts.csv"))

datanames <- c("C96_1","C96_2","S96_1","S96_2")

list <- llist(C96_1,C96_2,S96_1,S96_2)

# Make count matrix (keep gene_id as rownames and rawcounts)

fun <- function(i){
  i <- i %>% as.data.frame()
  rownames(i) <- i$ensembl_gene_id
  i <- i %>% dplyr::select(RawCounts, ensembl_gene_id)
}

list <- lapply(list,fun)
list2env(list, envir=.GlobalEnv)

head(C48_1)

# Change rawcounts to sample name


colnames(C96_1)[1] <- datanames[1]
colnames(C96_2)[1] <- datanames[2]
colnames(S96_1)[1] <- datanames[3]
colnames(S96_2)[1] <- datanames[4]

# Join the data tables to form the final matrix

counts <- list(C96_1,C96_2,S96_1,S96_2) %>% 
  purrr::reduce(dplyr::full_join,)

rownames(counts) <- counts$ensembl_gene_id
counts <- counts %>% dplyr::select(-ensembl_gene_id)
head(counts)

# Create the coldata for the high level data structure

coldata <- data.frame(
  celltype=c(rep("ctrl",2),rep("selected",2)),
  timeline = as.factor(rep(96,4)),
  replicate=rep(c(1,2),2))

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

# PCA preparation

rld <- rlogTransformation(dds)

plotPCA(rld,intgroup="celltype")  
plotPCA(rld,intgroup="timeline")
plotPCA(rld,intgroup="replicate") 

# Results can be extracted

res <- results(dds, name = "celltype_selected_vs_ctrl")

res_tbl <- as_tibble(res, rownames="ENSMUG")

# Save the tibbles and R objects

saveRDS(dds, file = "./data/221129_dds.rds")
saveRDS(res_tbl, file = "./data/221129_res_tbl.rds")
write_tsv(res_tbl, file = "./data/221129_res_tbl.tsv")

# Check the raw counts of the Slc7a5 channel in each sample

as_tibble(counts(dds["ENSMUSG00000040010"], normalize = TRUE),
          rownames = 'ENSMUG') %>%
  gather(sample, counts, -ENSMUG) %>%
  left_join(as_tibble(coldata, rownames = "sample")) %>%
  ggplot(aes(x = sample, y = counts, fill = celltype)) +
  geom_bar(stat = 'identity', color = "gray30", width = .75) +
  theme(axis.text.x = element_text(size = 11, angle = 1),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11))+
  ggtitle("Slc7a5 mRNA counts")

# Data annotation---------------------------------------------------------------

# Import annotation file

ensembl_to_geneName <- read_csv("./data/Biomart_annotations_mm10.csv")

# Add gene names and entrez id to your results

res_tbl <- res_tbl %>%
  left_join(ensembl_to_geneName) %>%
  arrange(padj)

write.csv(res_tbl,"./data_output/sel_vs_ctrl_T96/res_tbl.csv")

## Volcano plot-----------------------------------------------------------------

library(ggrepel)

target <- read_xlsx("./data/candidates.xlsx", sheet = 1)
target <- target$Gene
highlight <- signif %>% filter(gene%in%target)

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|padj<0.05&log2FoldChange<=-1,as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  geom_text_repel(max.overlaps = 50)+
  labs(title = "Selected vs control at T48")+
  theme_bw()

ggsave("./figures/sel_vs_ctrl_T48/volcanoplot", last_plot(), device = png, dpi= 500)


res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|padj<0.05&log2FoldChange<=-1,as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "lightcoral")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  theme_bw()+
  geom_point(data=highlight, 
             aes(x=log2FoldChange,y=-log10(padj)), 
             color='blue',
             size=0.75)+
  geom_text_repel(data = highlight, size = 4, segment.color = "blue",
                  max.overlaps = 20, min.segment.length = 0, color = "blue",
                  box.padding = 0.5)+
  theme(legend.title= element_blank())+
  labs(title = "Selected vs control at T48")

ggsave("./figures/sel_vs_ctrl_T48/volcanoplot_targets", last_plot(), device = png, dpi= 500)


## OVER REPRESENTATION ANALYSIS-------------------------------------------------

# GO and KEGG using ORA

res_tbl <- read_csv("./data_output/sel_vs_ctrl_T96/res_tbl.csv")

# Filter the significantly DE genes

diff <- res_tbl %>% 
  filter(abs(logFc)>=1&padj<0.05)

# ORA
# These terms are classed into three categories, called namespaces:

# Molecular Function (MF): molecular activities of gene products
# Cellular Component (CC): where gene products are active
# Biological Process (BP): pathways and larger processes made up of the 
# activities of multiple gene products

library("GO.db")
library("org.Mm.eg.db")
library("clusterProfiler")
library("enrichplot")

# GO

de_genes <- unique(diff$ENTREZID)
all_genes <- unique(res_tbl$ENTREZID)

go_ora <- enrichGO(gene = de_genes,
                   OrgDb = org.Mm.eg.db,
                   universe = all_genes,
                   ont = "ALL",
                   readable = TRUE) 

go_ora_tbl <- as_tibble(go_ora) # Reduce results by beong more stringent on thresholds


go_ora_cc <- go_ora%>% filter(ONTOLOGY=="CC")
go_ora_cc_tbl <- as_tibble(go_ora_cc)
go_ora_bp <- go_ora%>% filter(ONTOLOGY=="BP")
go_ora_bp_tbl <- as_tibble(go_ora_bp)
go_ora_mf <- go_ora%>% filter(ONTOLOGY=="MF")
go_ora_mf_tbl <- as_tibble(go_ora_mf)


barplot(go_ora, showCategory=40) + ggtitle("barplot for ORA")
ggsave("./figures/sel_vs_ctrl_T96/barplot go ORA all", plot = last_plot(), device = png, dpi = 400)
barplot(go_ora_cc, showCategory=40) + ggtitle("barplot for ORA cellular component")
ggsave("./figures/sel_vs_ctrl_T96/barplot go ORA cellular component", plot = last_plot(), device = png, dpi = 400)
barplot(go_ora_bp, showCategory=40) + ggtitle("barplot for ORA biological process")
ggsave("./figures/sel_vs_ctrl_T96/barplot go ORA biological process", plot = last_plot(), device = png, dpi = 400)
barplot(go_ora_mf, showCategory=40) + ggtitle("barplot for ORA molecular function")
ggsave("./figures/sel_vs_ctrl_T96/barplot go ORA molecular function", plot = last_plot(), device = png, dpi = 400)

dotplot(go_ora, showCategory=25) + ggtitle("dotplot for ORA")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go ORA all", plot = last_plot(), device = png, dpi = 400)
dotplot(go_ora_bp, showCategory=25) + ggtitle("dotplot for ORA bp")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go ORA bp", plot = last_plot(), device = png, dpi = 400)
dotplot(go_ora_cc, showCategory=25) + ggtitle("dotplot for ORA cc")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go ORA cc", plot = last_plot(), device = png, dpi = 400)
dotplot(go_ora_mf, showCategory=25) + ggtitle("dotplot for ORA mf")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go ORA mf", plot = last_plot(), device = png, dpi = 400)

# Datasets to query panther or gorilla or other
# Allows for external validation of our go_ora tibble

write.table(res_tbl$gene,"./data_output/sel_vs_ctrl_T96/backgrd_gene_names.txt", sep = "\t",
            row.names = FALSE, quote = F)
write.table(diff$gene, file = "./data_output/sel_vs_ctrl_T96/signif_gene_names.txt", sep = "\t",
            row.names = FALSE, quote = F)


# GSEA--------------------------------------------------------------------------

ordered_genes <- abs(res_tbl$log2fc) # could also use log2fc to order them
names(ordered_genes) <- res_tbl$ENTREZID
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

# GO

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Mm.eg.db,
                 scoreType = "pos",
                 ont          = "ALL",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

go_gsea_tbl <- as_tibble(go_gsea)

go_gsea

write.csv(go_gsea_tbl,"./data_output/sel_vs_ctrl_T96/go_gsea_tbl.csv")

dotplot(go_gsea, showCategory=30) + ggtitle("Dotplot for GSEA all")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA all", plot = last_plot(), device = png, dpi = 400)
dotplot(go_gsea %>% filter(ONTOLOGY=="BP"), showCategory=30)+ ggtitle("Dotplot for GSEA BP")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA BP", plot = last_plot(), device = png, dpi = 400)
dotplot(go_gsea %>% filter(ONTOLOGY=="CC"), showCategory=30)+ ggtitle("Dotplot for GSEA CC")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA CC", plot = last_plot(), device = png, dpi = 400)
dotplot(go_gsea %>% filter(ONTOLOGY=="MF"), showCategory=30)+ ggtitle("Dotplot for GSEA MF")
ggsave("./figures/sel_vs_ctrl_T96/dotplot go GSEA MF", plot = last_plot(), device = png, dpi = 400)

