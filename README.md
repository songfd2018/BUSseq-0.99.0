# Batch Effects Correction With Unknow Subtypes for scRNA-seq Data (BUSseq)

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for use](#instructions-for-use)
- [License](./LICENSE)
- [Citation](#citation)

# Overview
Single-cell RNA-sequencing (scRNA-seq) technologies enable the measurement of the transcriptome of individual cells, which provides unprecedented opportunities to discover cell types and understand cellular heterogeneity. Despite their widespread applications, single-cell RNA-sequencing (scRNA-seq) experiments are still plagued by batch effects and dropout events.

One of the major tasks of scRNA-seq experiments is to identify cell types for a population of cells. Therefore, the cell type of each individual cell is always unknown and is the target of inference. However, most existing methods for batch effects correction, such as Combat space and the surrogate variable analysis (SVA), are designed for bulk experiments and require knowledge of the subtype information, which corresponds to cell type information for scRNA-seq data, of each sample a priori.
  
Here, the R package `BUSseq` fits an interpretable Bayesian hierarchical model---the Batch Effects Correction with Unknown Subtypes for scRNA seq Data(BUSseq)---to correct batch effects in the presence of unknown cell types. BUSseq is able to simultaneously correct batch effects, clusters cell types, and takes care of the count data nature, the overdispersion, the dropout events, and the cell-specific sequencing depth of scRNA-seq data. After correcting the batch effects with BUSseq, the corrected value can be used for downstream analysis as if all cells were sequenced in a single batch. BUSseq can integrate the read count matrices measured from different platforms and allow cell types to be measured in some but not all of the batches as long as the experimental design fulfills the conditions listed in our [manuscript](https://www.biorxiv.org/content/10.1101/533372v1).

# Repo Contents

- [R](./R): `R` package code.
- [data](./data): example data for the demo.
- [inst/doc](./inst/doc): package pdf user's guide, and usage of the `BUSseq` package on a small example dataset.
- [man](./man): package manual for help in R session.
- [src](./src): `C++` source code.
- [tests](./tests): tests written using the `BUSseq` package.
- [vignettes](./vignettes): `R` vignettes for R session html help pages.

# System Requirements

## Hardware Requirements

The `BUSseq` package works on a standard personal computer (PC). The runtimes reported below were generated on a Windows 10 operating system by a PC with 16 GB RAM and 2 cores of 3.4 GHz.

## Software Requirements

### OS Requirements

The package supports *Linux*, *Mac* and *Windows* operating systems. It has been tested on the following systems:

Linux: Ubuntu 18.04

Mac OSX: Mac OS X 10.14 Mojave

Windows: Windows 10 Enterprise

### Software dependencies

Before installing the `BUSseq` package, users should have installed `R` with version 3.5.0 or higher. For Windows system, the users should also install [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/).

#### Installing R version 3.5.2 on Windows

Please download and install [R-3.5.2 for Windows](https://cran.r-project.org/bin/windows/base/) onto your computer. It will take a few minutes.

#### Package dependencies

Users should install the following packages prior to installing `BUSseq`, from an `R` session:

```
install.packages(c('devtools', 'bigmemory', 'gplots', 'knitr'))
```

which will install in about one minute.

#### Package Versions

The `BUSseq` package functions require the above packages with the following versions, respectively:

```
devtools: 2.0.1
bigmemory: 4.5.33
gplots: 3.0.1.1
knitr: 1.12
```

If you encounter any problem with installation, please drop us an [Issue](https://github.com/songfd2018/BUSseq/issues). 

# Installation Guide

From an `R` session, type:

```
require(devtools)
install_github("songfd2018/BUSseq", build_vignettes = TRUE) # install BUSseq
```

or

```
install.packages("/your/local/directory/BUSseq_0.99.6.tar.gz", repos = NULL, type = "source") # install BUSseq from zip file
```

It takes approximately 30 seconds to install directly from Github and it costs 10 seconds to install the compiled package from a local directory. 

# Demo

Please find the pdf file for the detailed instructions on how to use the package by running the following codes in the `R` session:

```
vignette("BUSseq_user_guide",package="BUSseq")  # view the vignettes
```

For a given number of cell types *K*, it takes about 6 minutes to run on the demo dataset. To select the optimal *K*, we need to compare the BIC values for different *K*s. In the vignettes, we enumerate *K* from 2 to 6. As a result, we need about 5 * 6 = 30 minutes. When we have a multi-core computer or a cluster, we can run BUSseq with different *K*s in parallel. In that case, with 5 cores, we can finish all the computation within 10 minutes.

# Instructions for use
Please follow the steps in the [BUSseq implementation](https://github.com/songfd2018/BUSseq_implementation) repository to reproduce all the results and figures of simulation and real data analysis.

# Citation
BUSseq manuscript is available from [bioRxiv](https://www.biorxiv.org/content/10.1101/533372v1). If you use BUSseq for your work, please cite our paper.

		@article {Song533372,
			author = {Song, Fangda and Chan, Ga Ming and Wei, Yingying},
			title = {Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effects Correction},
			elocation-id = {533372},
			year = {2019},
			doi = {10.1101/533372},
			publisher = {Cold Spring Harbor Laboratory},
			URL = {https://www.biorxiv.org/content/early/2019/01/29/533372},
			eprint = {https://www.biorxiv.org/content/early/2019/01/29/533372.full.pdf},
			journal = {bioRxiv}
		}
