# MICA: A multi-omics method to predict gene regulatory networks in early human embryos

Gregorio Alanis-Lobato<sup>1</sup>, Thomas E. Bartlett<sup>2</sup>, Qiulin Huang<sup>1,3</sup>, Claire Simon<sup>1</sup>, Afshan McCarthy<sup>1</sup>, Kay Elder<sup>4</sup>, Phil Snell<sup>4</sup>, Leila Christie<sup>4</sup> & Kathy Niakan<sup>1,3</sup>

1. Human Embryo and Stem Cell Laboratory, The Francis Crick Institute, London, UK
2. Department of Statistical Science, University College, London WC1E 7HB, UK
3. The Centre for Trophoblast Research, Department of Physiology, Development and Neuroscience, University of Cambridge, Cambridge CB2 3EG, UK.
4. Bourn Hall Clinic, Bourn, Cambridge, CB23 2TN.

## Introduction

Under construction...

## Repo organisation

This repository is organised as follows:

- `01_data_preprocessing`: Scripts used to preprocess raw scRNA-seq and liATAC-seq data (FASTQ files).
- `02_data_processing`: Scripts for QC, filtering and normalisation of scRNA-seq and liATAC-seq data.
- `03_grn_inference`: Implementation of the different methods for GRN inference.
- `04_grn_evaluation`: Scripts for statistical and biological evaluation of the inferred GRNs.
- `05_figures`: Scripts to construct the main and supplementary figures from the manuscript.
