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

One of the major tasks of scRNA-seq experiments is to identify cell types for a population of cells. Therefore, the cell type of each individual cell is always unknown and is the target of inference. However, most existing methods for batch effects correction, such as Combat and the surrogate variable analysis (SVA), are designed for bulk experiments and require knowledge of the subtype information, which corresponds to cell type information for scRNA-seq data, of each sample a priori.
  
Here, the R package `BUSseq` fits an interpretable Bayesian hierarchical model---the Batch Effects Correction with Unknown Subtypes for scRNA seq Data (BUSseq)---to correct batch effects in the presence of unknown cell types. BUSseq is able to simultaneously correct batch effects, clusters cell types, and takes care of the count data nature, the overdispersion, the dropout events, and the cell-specific sequencing depth of scRNA-seq data. After correcting the batch effects with BUSseq, the corrected value can be used for downstream analysis as if all cells were sequenced in a single batch. BUSseq can integrate read count matrices obtained from different scRNA-seq platforms and allow cell types to be measured in some but not all of the batches as long as the experimental design fulfills the conditions listed in our [manuscript](https://www.biorxiv.org/content/10.1101/533372v1).

# Repo Contents

- [R](./R): `R` code.
- [data](./data): the example data for the demo.
- [inst/doc](./inst/doc): compiled user's guide and illustration of the applications of `BUSseq` package to the demo dataset.
- [man](./man): help manual.
- [src](./src): `C++` source code.
- [tests](./tests): sample code for the demo dataset.
- [vignettes](./vignettes): source code for the user's guide.

# System Requirements

## Hardware Requirements

The `BUSseq` package works on a standard personal computer (PC). The runtimes reported below were generated on an Ubuntu 18.04 operating system by a PC desktop with 8 GB RAM and 8 cores of 2.6 GHz.

## Software Requirements

### OS Requirements

The package supports *Linux*, *Mac* and *Windows* operating systems. It has been tested on the following systems:

Linux: Ubuntu 18.04

Mac OSX: Mac OS X 10.14 Mojave

Windows: Windows 10 Enterprise

### Software dependencies

Before installing the `BUSseq` package, users should have installed `R` with version 3.5.0 or higher. For Windows system, the users should also install [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/).

#### Installing R version 3.5.2 on Ubuntu 18.04

To install the latest stable version of R on Ubuntu 18.04, please follow these steps:

```
sudo apt install apt-transport-https software-properties-common
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
sudo apt update
sudo apt install r-base
R –version
```

The installation will take a few minutes.

#### Package dependencies

Users should install the following packages prior to installing `BUSseq`, from an `R` session:

```
install.packages(c('devtools', 'bigmemory', 'gplots', 'knitr'))
```

The installation will take about one minute.

#### Package Versions

The `BUSseq` package depends on the above packages with the following versions, respectively:

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

Please check the user's guide for the detailed instructions on how to use the package by running the following code in the `R` session:

```
vignette("BUSseq_user_guide",package="BUSseq")  # view the vignettes
```

For a given number of cell types *K*, it takes about 3 minutes to run the analysis for the demo dataset. To select the optimal *K*, we need to compare the BIC values for different *K*s. In the vignettes, we enumerate *K* from 2 to 6. As a result, we need about 5 * 3 = 15 minutes. When we have a multi-core computer or a cluster, we can further run BUSseq with different *K*s in parallel. For the simulation study in our manuscript which consists of 10,000 genes and 600 cells, it takes 4.25 hours to run 3,000 iterations of the MCMC algorithm for a given K. 


# Instructions for use
Please follow the steps in the [BUSseq implementation](https://github.com/songfd2018/BUSseq_implementation) repository to reproduce all the results and figures in our [manuscript](https://www.biorxiv.org/content/10.1101/533372v1).

BUSseq stores the MCMC samples onto the hard disk to reduce memory usage. For the  simulation dataset, mouse hematopoietic and human pancreatic studies in our manuscript, each takes 6~8 GB to save the MCMC samples. Please ensure that there is enough space on your hard disk before running the codes for reproduction.

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