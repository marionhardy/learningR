library(tidyverse)
library(GEOquery)

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


counts

counts[-1,]
coln <- counts[1,]

colnames(counts)<- c("PBS_1","PBS_2","PBS_3",
                     "231_1","231_2","231_3",
                     "NME1_1","NME1_2","NME1_3")

head(counts)

counts <- counts[-1,]

class(counts$PBS_1)
counts <-  type.convert(counts, as.is = TRUE) 


coldata <- data.frame(
  celltype = rep("MCF7",9),
  condition = c(rep("PBS",3),rep("231",3),rep("NME1",3)),
  replicates = as.factor(c(rep(1:3,3))))

rownames(coldata)<- c("PBS_1","PBS_2","PBS_3",
                      "231_1","231_2","231_3",
                      "NME1_1","NME1_2","NME1_3")


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~ condition ) 

dds <- DESeq(dds)

sizeFactors(dds) 

plotDispEsts(dds) 

resultsNames(dds) 

dds$condition <- relevel(dds$condition, ref="PBS")
dds <- DESeq(dds)
resultsNames(dds) 

res <- results(dds, name = "condition_NME1_vs_PBS")
res_tbl <- as_tibble(res, rownames="Gene_names")

# PCA to check experimental design

rld <- rlogTransformation(dds)

plotPCA(rld,intgroup="condition")  
plotPCA(rld,intgroup="replicates") 

# Save the tibbles etc

saveRDS(dds, file = "./data/GSE189685_dds.rds")

# Volcano plot

library(ggrepel)

# target <- read_xlsx("./data/candidates.xlsx", sheet = 1)
# target <- target$Gene
# highlight <- signif %>% filter(gene%in%target)

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|padj<0.05&log2FoldChange<=-1,
                            as.character(Gene_names),''))) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  geom_text_repel(max.overlaps = 50)+
  labs(title = "NME1 compared to PBS treatment")+
  theme_bw()

ggsave("./figures/volcano_nme1vsPBS.png", last_plot(), device = png, dpi= 500,
       width = 12, height = 8)


# Volcano plot of lfc >-0.5 <0.5 , color blue
# Lfc 0.264 means ? 1,2x change in gene expression

res_tbl %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 0.264, 
             label = ifelse(padj<0.05&log2FoldChange>=0.264|padj<0.05&log2FoldChange<=-0.264,
                            as.character(Gene_names),''))) +
  scale_colour_manual(values = c("gray", "blue")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 0.264) +
  geom_vline(xintercept = -0.264)+
  geom_text_repel(max.overlaps = 30, size = 3)+
  labs(title = "NME1 compared to PBS treatment")+
  theme_bw()

ggsave("./figures/volcano_nme1vsPBS_415.png", last_plot(), device = png, dpi= 500,
       width = 12, height = 8)

# How many genes are DEG under our conditions

res_tbl %>% 
  filter(abs(log2FoldChange)>=0.264, padj<0.05) 

# 438 genes

## OVER REPRESENTATION ANALYSIS-------------------------------------------------

# If we remember well, the input to do ORA are ENTREZID. But the function has a 
# specific argument called keytype which might have other options such as
# gene names. This would allow us to not have to annotate our results with
# ENTREZ ids

# GO and KEGG using ORA

# Filter the significantly DE genes

diff_1 <- res_tbl %>% 
  filter(abs(log2FoldChange)>=1&padj<0.05)

diff_0264 <- res_tbl %>% 
  filter(abs(log2FoldChange)>=0.264&padj<0.05)

# ORA

library("GO.db")
# library("org.Mm.eg.db") ! we're in a human organism
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")

# GO

res_tbl <- drop_na(res_tbl)

de_genes <- unique(diff_1$Gene_names) # only 8 differentialy expressed genes
all_genes <- unique(res_tbl$Gene_names)

keytypes(org.Hs.eg.db)

go_ora <- enrichGO(gene = de_genes,
                   keyType = "GENENAME",
                   OrgDb = org.Hs.eg.db,
                   universe = all_genes,
                   ont = "ALL",
                   readable = TRUE) # it's empty

de_genes <- unique(diff_0264$Gene_names) # 438 DEGs
all_genes <- unique(res_tbl$Gene_names)

go_ora <- enrichGO(gene = de_genes,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   universe = all_genes,
                   pAdjustMethod = "BH",
                   ont = "ALL",
                   readable = TRUE) 
go_ora

barplot(go_ora, showCategory=40) + ggtitle("barplot for ORA")
ggsave("./figures/ORA_ALL.png", plot = last_plot(), device = png, dpi = 300,
       width = 12, height = 18)

barplot(go_ora %>% filter(ONTOLOGY=="MF"), showCategory=20)

# GSEA--------------------------------------------------------------------------

ordered_genes <- abs(res_tbl$log2FoldChange)
names(ordered_genes) <- res_tbl$Gene_names
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

# GO

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Hs.eg.db,
                 scoreType = "pos",
                 keyType = "SYMBOL",
                 ont          = "ALL",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

dotplot(go_gsea, showCategory=30) + ggtitle("Dotplot for GSEA all")
ggsave( "./figures/GSEA_ALL.png", plot = last_plot(), device = "png", dpi = 300,
       width = 8, height = 7)

dotplot(go_gsea %>% filter(ONTOLOGY=="MF"), showCategory=30)+ ggtitle("Dotplot for GSEA BP")
ggsave("./figures/GSEA_MF.png", plot = last_plot(), device = png, dpi = 300,
       width = 6, height = 9)

dotplot(go_gsea, split = "ONTOLOGY",showCategory=10, x="NES") + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all split")
ggsave("./figures/GSEA_ALL_split.png", plot = last_plot(), device = png, dpi = 300,
       width = 8, height = 12)

go_gsea_tbl <- as.tibble(go_gsea)

# What is interesting?

# There's a lot of cell cycle activity in cells that were treated with the H1 exo
# There's response to inflammation 
# modifications in channel activity GO:0015267 
# signaling receptor GO:0038023

# Check for specific pathways using KEGG and wikipathways

library(msigdbr)
print(msigdbr_collections(), n=23)

# Database containing annotated gene sets that can be used for pathway or gene set analyses 
# They have 9 collections : Hallmark, C1-C8
# https://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#H

hsa_reactome_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "C2",
  subcategory = "CP:REACTOME") # for reactome collection

hsa_kegg_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "C2",
  subcategory = "CP:KEGG") # for KEGG collection

hsa_wiki_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "C2",
  subcategory = "CP:WIKIPATHWAYS") # for Wikipathways collection

set.seed(69)


gsea_results_react <- GSEA(
  geneList = ordered_genes, 
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    hsa_reactome_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
) # There's IFN beta signaling

dotplot(gsea_results_react, x = "NES", showCategory = 30)+ ggtitle("GSEA reactome LFC")
ggsave("./figures/dotplot_GSEA_reactome.png", plot = last_plot(), device = png, dpi = 400,
       width = 10, height = 8)

gsea_results_kegg <- GSEA(
  geneList = ordered_genes,
  pvalueCutoff = 0.05, 
  eps = 0,
  seed = TRUE, 
  pAdjustMethod = "BH",
  scoreType = "pos",
  TERM2GENE = dplyr::select(
    hsa_kegg_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_kegg, x = "NES", showCategory = 30)+ ggtitle("GSEA reactome LFC")
ggsave("./figures/dotplot_GSEA_kegg.png", plot = last_plot(), device = png, dpi = 400,
       width = 10, height = 8)

# It's still possible to use pathview to extract the significant paths

library("pathview")

keggresids <- c("04350","04060")

foldchanges <- diff_0264$log2FoldChange
names(foldchanges) <- diff_0264$Gene_names
head(foldchanges)
table(is.na(foldchanges))

tmp <- sapply(keggresids, function(pid) pathview(gene.data = foldchanges,
                                                 gene.idtype = "SYMBOL",
                                                 pathway.id = pid,
                                                 species = "hsa",
                                                 kegg.dir="./figures/",
                                                 out.suffix= "_Colored",
                                                 kegg.native = TRUE,
                                                 map.null = FALSE))


gsea_results_wiki <- GSEA(
  geneList = ordered_genes, 
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hsa_wiki_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_wiki, x = "NES", showCategory = 30)+ ggtitle("GSEA Wikipathway LFC")

ggsave("./figures/dotplot_GSEA_Wikipathway.png", plot = last_plot(), device = png, dpi = 400)




## NB their processed data file is an intermediate step in the processing of DEseq2
# i believe. It cannot be used as input to DESeq2 since it comprises integers









