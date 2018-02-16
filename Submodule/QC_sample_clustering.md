## Learning Objectives

* Understand RNA-seq data normalization
* Create a `DESeqDataset` object
* Understand the priciples of PCA and Hierarchical Clustering
* Perform PCA and HC on our data

## Normalization

Normalization of reads is a necessary step for all visualization procedures, including PCA and hierarchical clustering. **The DE analysis in both DESeq2 and edgeR takes raw counts, and performs their own normalization as part of the DE analysis, normalization is for visualization only (or other downstream analyses)**.




