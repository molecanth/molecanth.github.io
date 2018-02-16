## Learning Objectives

* Set up a DE project in R 
* Understand the use of negative binomial distribution to model RNA-seq count data
* Perform DE Analysis

## Set up a DE analysis in R

Before we begin we will nedd to set up a new project for this analysis in RStudio: 

1. Go to the `File` menu and select `New Project`.
2. In the `New Project` window, choose `New Directory`. Then, choose `Empty Project`. Name your new directory `DEanalysis` and then "Create the project as subdirectory of:" the Desktop (or location of your choice).
3. The new project should automatically open in RStudio. 

To check whether or not you are in the correct working directory, use `getwd()`. The path `Desktop/DEanalysis` should be returned to you in the console. Within your working directory use the `New folder` button in the bottom right panel to create three new directories: `data`, `meta` and `results`. An important key to a good analysis is establishing project organization before you start. 

Go to the `File` menu and select `New File`, then select `R Script`. This should open up a script editor in the top left hand corner. This is where we will be typing and saving all commands required for this analysis. In the script editor type in header lines:

```
## Gene-level differential expression analysis using DESeq2
```

Now save the file as `de_script.R`. When finished your working directory should now look similar to this:

![setup](../img/settingup.png)

Finally, we need to grab the files that we will be working with for the analysis. The NHPRTR has provided a massive count matrix of every species/tissue combination, of which I have extracted a subset of data for skeletal muscle from all species, and generated an associated metadata file. Right click on the links below, and choose the "Save link as ..." option to download:

* Save the [full counts matrix](https://github.com/molecanth/molecanth.github.io/blob/master/module_data/primate_skeletalmuscle.txt) file in the `data` directory.
* Save the [full metadata table](https://github.com/molecanth/molecanth.github.io/blob/master/module_data/primate_meta.txt) file in the `meta` directory.

## Loading libraries

For this analysis we will be using several R packages, some which have been installed from CRAN and others from Bioconductor. To use these packages (and the functions contained within them), we need to **load the libraries.** Add the following to your script and don't forget to comment liberally!

```r
## Setup
### Bioconductor and CRAN libraries used
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
```

## Loading data

To load the data into our current environment, we will be using the `read.table` function. We need to provide the path to each file and also specify arguments to let R know that we have a header (`header = T`) and the first column is our row names (`row.names =1`). By default the function expects tab-delimited files, which is what we have.

```r
## Load in data
data <- read.table("data/Mov10_full_counts.txt", header=T, row.names=1) 

meta <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)
```

Use `class()` to inspect our data and make sure we are working with data frames:

```r
### Check classes of the data we just brought in
class(meta)
class(data)
```

## Viewing data

Make sure your datasets contain the expected samples / information before proceeding to perfom any type of analysis. 

```r
View(meta)
View(data)
```
## DE analysis overview

So what does this count data actually represent? The count data used for differential expression analysis represents the number of sequence reads that originated from a particular gene. The higher the number of counts, the more reads associated with that gene, and the assumption that there was a higher level of expression of that gene in the sample. 

<img src="../img/deseq_counts_overview.png" width="600">

With DE analysis, we are looking for genes that change in expression between two or more groups (defined in the metadata)with respect to both the biological and technical variation between conditions, for example:
- treatment vs. control

**Why does it not work to identify differentially expressed gene by ranking the genes by how different they are between the two groups based only on fold change values?**

All RNA-seq datasets have variation between conditions that is not a result the conditions themselves. The goal of DE analysis to indentify genes that are differentially expressed between conditions while controlling for variation introduced by variables that are not of interest.

Even though the mean expression levels between sample groups may appear to be quite different, it is possible that the difference is not actually significant when technical variation is not controlled for. This is illustrated for 'GeneA' expression between 'untreated' and 'treated' groups in the figure below. The mean expression level of geneA for the 'treated' group is twice as large as for the 'untreated' group, but the variation between replicates indicates that this may not be a significant difference. While it's possible that your treatment is driving the variation in transcription, this situation or more often the result of some technical variable, such as a batch effect. **We need to take into account the variation in the data (and where it might be coming from) when determining whether genes are differentially expressed.**

<img src="../img/de_norm_counts_var.png" width="400">





