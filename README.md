# learningR

Learning R basics for an intro to RNAseq data processing and analysis.

Objectives: Teaching people who are not familiar with any coding interface to make plots, use DESeq2 and do enrichment analyses using MSigDb on basic RNAseq data + help them read error messages and find resources to troubleshoot. Course given to the BVDE lab during 2022-2023.

Caveats: Not enough time to make people truly familiar with R and be able to write code from scratch. Additionally, it is advised to consult multiple different resources about RNAseq data analysis and statistics in parallel to this course. Some are listed at the end of the ch08 powerpoint.

### Ch01 Set up

- Setting up R, RStudio.
- Installing base packages and using biocmanager.

### Ch02 Introduction to R

- Good practices
- Data types and structures
- Practice

### Ch03 Data manipulation

- Modifying data structures
- Introduction to the tidyverse package: usinf filter(), mutate(), group_by() and summarize()
- Practice

### Ch04 Data manipulation pt 2

- Practice from Ch03 but using a dataset from gene expression omnibus

### Ch05 Data visualization

- Using ggplot2()
- Using plots for data quality control
- Using appropriate plots and visuals to answer questions about data

### Ch06 Q&A session + intro to high level data structures

- Violinplot and adding conditional highlights
- Types of high level structures and their use

### Ch07 DESeq2 (accompanied by the Ch08 powerpoint!)

- DESeq2 object creation
- Writing and using a loop function on a dataframe
- Quality check metrics : size factors, dispersion estimation and pca

### Ch08 RNAseq: processing and analysis pipeline

Powerpoint presentation covering:
- Data processing (FASTQ -> BAM -> counts) using trim_galore, RSubread and featureCounts
- Data analysis for differential expression using DESeq2 (+ which package to choose and a statistics reminder)
- ORA vs GSEA
- Querying MSidDb for enrichment analyses

### Ch09 RNAseq: from Gene Expression Omnibus to a volcano plot

- Practice what was learned in ch07+08 on data from Gene Expression Omnibus (+ how to dowNload it)
- Questions are in the "Learning RNAseq" document

NB: As of 12/12/22, the article from which originated the data is up for review so the data is not available through GEO. To receive the dataset used for the practice, contact me. (the article is planned to be published on 12/12/23).


