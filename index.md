## UO Molecular Anthropology Group Genomics Training Module

This module is an end-to-end differential expression analysis (DGE) workflow which can be applied to a variety of biological questions using RNA-seq data. The module requires a working knowledge of R and will utilize the R integrated environment RStudio. 

Go to [Submodule 1](Submodules

### Installation Requirements

Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) 
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
 
Install the required R packages by running the following code in RStudio:

```r
source("http://bioconductor.org/biocLite.R") 
biocLite(c("RColorBrewer", "pheatmap", "gProfileR", "DESeq2", "clusterProfiler", 
           "DOSE", "org.Hs.eg.db", "pathview", "treemap", "purrr", "SPIA", "DEGreport"))
```

Load the libraries to make sure the packages installed properly:

```r
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(SPIA)
library(purrr)
library(gProfileR)
library(treemap)
```

Begin with [DGE Overview](https://github.com/molecanth/molecanth.github.io/blob/master/DGE_overview.md)

You can use the [editor on GitHub](https://github.com/molecanth/molecanth.github.io/edit/master/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/molecanth/molecanth.github.io/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
