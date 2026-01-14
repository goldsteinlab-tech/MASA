# MASA Library

[![DOI](https://doi.org/)
[![Conda Version](https:)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


**MASA** MASA – a computational tool to predict transcription factor activity from ATAC-seq
MASA (Motif-Associated Signs of Activity) is a computational framework that infers differential TF activity directly from ATAC-seq data in an unbiased and multilayered manner. 

## Abstract

Transcription factors (TFs) are central mediators of gene regulation, interpreting upstream signals into context-specific transcriptional programs. 
Chromatin accessibility profiling by ATAC-seq provides a powerful, genome-wide readout of regulatory activity, yet robust inference of TF activity from ATAC-seq remains challenging due to motif redundancy,
limited comparability across conditions, and reliance on a single analytical approach. Here, we present MASA (Motif-Associated Signs of Activity),
a computational framework that infers differential TF activity directly from ATAC-seq data in an unbiased and multilayered manner. 
MASA quantifies three complementary indicators of TF activity for each motif: footprint depth, flanking chromatin accessibility, and motif enrichment, and integrates them into an intuitive plot and comprehensive result tables.
To reduce redundancy and improve interpretability, MASA operates on archetypal consensus motif clusters rather than individual, highly similar motifs.
MASA requires only ATAC-seq data, does not depend on prior knowledge of relevant TFs, and does not necessitate focusing exclusively on differentially accessible regions.
We validate MASA in both mouse and human datasets across multiple biological contexts, including nutritional transitions, genetic perturbation, and human disease states.
In each case, MASA accurately predicts TFs with altered activity, which we independently confirm using ChIP-seq data.
In each case, MASA accurately predicts motifs with altered activity and nominates candidate TFs, which we confirm using representative TF occupancy.
MASA provides a robust solution for extracting biologically meaningful TF activity signatures from chromatin accessibility data and is well suited for diverse experimental designs.

## Installation

## Overview

**masa** is a publicly released R package.
The code in this repository corresponds to the final version described in the associated publication(s).

### Using R
The current official release of **masa** is available on GitHub.

You can install it in R with:

```r
# install.packages("remotes") # run this once if you don't have it
remotes::install_github("goldsteinlab-tech/masa", ref = "main")
```

### Using conda 

```bash
conda install -c conda-forge masa
```
Using pip
```bash
pip install masa
```

### Requirements
R >=4.3

DEseq2 R library

HOMER v5

#### Using conda:
python >=3.14


Tested on Linux 

### Quick Start
See the vignette

### Documentation
Full documentation, including tutorials and API reference, is available at:
[https://goldsteinlab-tech.github.io/MASA/docs/masa_vignette.html]

The documentation covers:

Installation and configuration

End-to-end tutorials

Reproducible analysis workflows

### Workflow reproducibility
The data to reproduce the vignette workflow are avalaible here:
ZELANDO

The motifDB is available here:
mm10:
hg19:

The nuccode files are avialable here:
mm10:
hg19:

The complete results from the workflow is available here:
DOI

### Citation
If you use MASA in scientific work, please cite:

text
@software{masa_2026,
  author    = {Leslie Cohen and Ido Goldstein},
  title     = {MASA – a computational tool to predict transcription factor activity from ATAC-seq},
  year      = {2026},
  publisher = {xxx},
  doi       = {xxx},
  url       = {xxx}
}
Zenodo record: https://doi.org/10.5281/zenodo.XXXXXXX


### License
This project is licensed under the GNU License – see the LICENSE file for details.

### Acknowledgments
This work was supported by the European Research Council to IG (ERC-StG #947907).
We thank Noga Korenfeld for her contributions and feedback.

### Contact
Maintainer: ido.goldstein@mail.huji.ac.il

Issues and bugs: https://github.com/goldsteinlab-tech/MASA/issues

