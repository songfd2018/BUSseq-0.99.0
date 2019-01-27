# BUSseq R package
The BUSseq R package implements the BUSseq model to adjust single-cell RNA-sequencing data for batch effects when there are unknown cell types.

## Introduction
Single-cell RNA-sequencing (scRNA-seq) technologies enable the measurement of the transcriptome of individual cells, which provides unprecedented opportunities to discover cell types and understand cellular heterogeneity [1]. Despite their widespread applications, single-cell RNA-sequencing (scRNA-seq) experiments are still plagued by batch effects and dropout events.

One of the major tasks of scRNA-seq experiments is to identify cell types for a population of cells [1]. Therefore, the cell type of each individual cell is always unknown and is the target of inference. However, most existing methods for batch effects correction, such as Combat space [2] and the surrogate variable analysis (SVA) [3, 4], are designed for bulk experiments and require knowledge of the subtype information, which corresponds to cell type information for scRNA-seq data, of each sample a priori.
  
Here, the R package \Rpackage{BUSseq} fits an interpretable Bayesian hierarchical model---the Batch Effects Correction with Unknown Subtypes for scRNA seq Data(BUSseq)---to correct batch effects in the presence of unknown cell types [5]. BUSseq is able to simultaneously correct batch effects, clusters cell types, and takes care of the count data nature, the overdispersion, the dropout events, and the cell-specific sequencing depth of scRNA-seq data. After correcting the batch effects with BUSseq, the corrected value can be used for downstream analysis as if all cells were sequenced in a single batch. BUSseq can integrate the read count matrices measured from different platforms and allow cell types to be measured in some but not all of the batches as long as the experimental design fulfills the conditions listed in [5].

## Installation
You can use the following command to install it.  

```
library(devtools)
install_github("songfd2018/BUSseq")
```

## User's Guide
Please refer to the vignetee for detailed function instructions using

```
browseVignettes("BUScorrect")
```

## References
1. Bacher, Rhonda, and Christina Kendziorski. "Design and computational analysis of single-cell RNA-sequencing experiments." *Genome Biology* 17, no. 1 (2016): 63.
2. Johnson, W. Evan, Cheng Li, and Ariel Rabinovic. "Adjusting batch effects in microarray expression data using empirical Bayes methods." *Biostatistics* 8.1 (2007): 118-127.
3. Leek, Jeffrey T., and John D. Storey. "Capturing heterogeneity in gene expression studies by surrogate variable analysis." *PLoS Genetics* 3.9 (2007): e161.
4. Leek, Jeffrey T. "Svaseq: removing batch effects and other unwanted noise from sequencing data." *Nucleic Acids Research* 42, no. 21 (2014): e161-e161.
5. Fangda Song, Ga Ming Chan and Yingying Wei. Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effects Correction, *Manuscript*.
