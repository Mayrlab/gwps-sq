[![DOI](https://zenodo.org/badge/532343780.svg)](https://zenodo.org/doi/10.5281/zenodo.10895729)

## Overview

This repository provides initial analysis of [K562 and RPE1 Perturb-seq data](https://gwps.wi.mit.edu) 
after reprocessing with the `scUTRquant` pipeline to resolve 3'UTR isoforms rather than gene counts. 
It includes both some preprocessing code used to created the metadata tables used in the scUTRquant 
pipeline, and the initial stages of filtering, summarization to pseudobulk, clustering of the KD6 
(K562 6-day essential screen) data, differential gene and isoform expression testing, and a replication
analysis using the RD7 (RPE1 7-day essential screen) data.

Additionally analyses downstream of this are included in [the `scUTRquant-figures` repository](https://github.com/Mayrlab/scUTRquant-figures).

The accompanying manuscript is openly available at:

> Fansler, M.M., Mitschka, S. & Mayr, C. Quantifying 3′UTR length from scRNA-seq data reveals changes independent of gene expression. *Nat Commun* **15**, 4050 (2024). [https://doi.org/10.1038/s41467-024-48254-9](https://doi.org/10.1038/s41467-024-48254-9)


## Organization
The folders in the repository have the following purposes:

- `analysis` - primary source code and rendered HTMLs of R Markdown
- `envs` - Conda environment YAML files for recreating the execution environment
- `img` - *output* images
- `metadata` - *output* tables used in pipeline
- `scripts` - additional scripts for data format conversions
- `tbl` - *output* result tables

All code is expected to be executed with this repository as the present working
directory. If opening as an R Project in RStudio, make sure to set the Project 
folder as the working directory.

### Source Code
The primary source code is found in the `analysis` folder. 
Files are numbered in the original order of execution, though the order does not 
imply strict necessity (most analyses here can be independently executed).

The `analysis/preprocessing/*_scutrquant_inputs.Rmd` files were used to populate the
`metadata` tables used by the scUTRquant pipeline. These should not be necessary to rerun
if starting scUTRquant outputs.

### Execution Environments
The R instances used to execute the files was captured both in the rendered RMDs themselves
(see **Runtime Details** section in HTMLs) and provided as YAML files in the `envs` folder.

To recreate on arbitrary platforms (Linux or MacOS), we recommend using 
[Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html#)
and the minimal YAML (`*.min.yaml`):

```bash
micromamba create -n bioc_3_16 -f envs/bioc_3_16.min.yaml
micromamba activate bioc_3_16
```

A fully-solved environment capture is also provided (`*.full.yaml`). This is only 
expected to recreate on the **osx-64** platform and is primarly intended for *exact* 
replication and a statement of record.
