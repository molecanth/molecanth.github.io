## Learning Objectives 

* Introduce DGE and experimental design
* Set up an RNA-seq project in RStudio 
* Understand RNA-seq and DGE analysis workflow
* Understand the use of negative binomial distribution to model RNA-seq count data

# Differential gene expression (DGE) analysis overview 

The goal of RNA-seq is often to perform differential expression testing to determine which genes are expressed at different levels between conditions. These genes can offer biological insight into the processes affected by the condition(s) of interest. 

To determine the expression levels of genes, our RNA-seq workflow followed the steps detailed in the image below. All steps were performed on the command line (Linux/Unix) through the generation of the read counts per gene. The differential expression analysis and any downstream functional analysis are generally performed in R using R packages specifically designed for the complex statistical analyses required to determine whether genes are differentially expressed.


<img src="../img/RNAseqWorkflow.png" width="600">

In the next few lessons, we will walk you through an **end-to-end gene-level RNA-seq differential expression workflow** using various R packages. We will start with the count matrix, perform exploratory data analysis for quality assessment and to explore the relationship between samples, perform differential expression analysis, and visually explore the results prior to performing downstream functional analysis.
