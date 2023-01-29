# MICA: A multi-omics method to predict gene regulatory networks in early human embryos

Gregorio Alanis-Lobato<sup>1</sup>, Thomas E. Bartlett<sup>2</sup>, Qiulin Huang<sup>1,3</sup>, Claire Simon<sup>1</sup>, Afshan McCarthy<sup>1</sup>, Kay Elder<sup>4</sup>, Phil Snell<sup>4</sup>, Leila Christie<sup>4</sup> & Kathy Niakan<sup>1,3</sup>

1. Human Embryo and Stem Cell Laboratory, The Francis Crick Institute, London, UK
2. Department of Statistical Science, University College, London WC1E 7HB, UK
3. The Centre for Trophoblast Research, Department of Physiology, Development and Neuroscience, University of Cambridge, Cambridge CB2 3EG, UK.
4. Bourn Hall Clinic, Bourn, Cambridge, CB23 2TN.

## Background

Recent advances in single cell -omics have been transformative in the characterisation of challenging-to-study biological contexts, including when the source material is extremely precious, such as in the early human embryo. Single cell datasets bring technical challenges to infer transcription factor-gene regulatory interactions, such as low read-depth leading to zero inflated data. Here we have systematically assessed the application of four different machine learning linear or non-linear gene regulatory network prediction strategies to single cell simulated and human embryo transcriptome datasets. We have also compared how the method of gene expression normalisation impacts on regulatory network predictions. Integrating chromatin accessibility datasets together with transcript expression datasets improved the reproducibility of the predicted gene regulatory networks. We found that the application of a non-linear network prediction method based on mutual information (MI) to single cell transcriptome datasets refined with chromatin accessibility (CA) data (called MICA), exhibited higher reproducibility compared to the alternative network interference methods tested. Moreover, MICA was used to make predictions about GRNs in the preimplantation human embryo, which were supported by previous findings in other developmental and stem cell contexts. Based on the gene regulatory networks predicted by MICA, we discovered co-localisation of the AP-1 transcription factor subunit proto-oncogene JUND and the TFAP2C transcription factor AP-2 in human preimplantation embryos. Overall, our comparative analysis of gene regulatory network prediction methods defines a pipeline that can be implemented on single-cell multi-omics datasets to infer interactions between transcription factor expression and target gene regulation that can then be functionally tested with laboratory interventions.

## Repo organisation

This repository is organised as follows:

- `01_data_preprocessing`: Scripts used to preprocess raw scRNA-seq and liATAC-seq data (FASTQ files).
- `02_data_processing`: Scripts for QC, filtering and normalisation of scRNA-seq and liATAC-seq data.
- `03_grn_inference`: Implementation of the different methods for GRN inference.
- `04_grn_evaluation`: Scripts for statistical and biological evaluation of the inferred GRNs.
- `05_figures`: Scripts to construct the main and supplementary figures from the manuscript.
