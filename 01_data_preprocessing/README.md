# Scripts used to preprocess raw scRNA-seq and liATAC-seq data (FASTQ files)

The following scripts and files are provided:

- `dw_fastq.R`: Provides functions to download FASTQ files from an FTP server given a CSV file with sample metadata and FTP addresses.
- `sb_dw_fastq`: Example `sbatch` script that calls `dw_fastq.R` to download FASTQ files in a computer cluster with the SLURM queueing system.
- `atacseq`
    - `design_hs_atac.csv`: CSV file with the list of liATAC-seq samples used in this study.
    - `nf_atac.sh`: Nextflow call used to preprocess the liATAC-seq data with the nf-core/atacseq v1.1.0 pipeline.
    - `regulatory_genomics_toolbox.yml`: YAML file with the instructions to create the conda environment used for TF footprinting and TF motif enrichment analysis.
    - `diff_footprints.sh`: Example script to use the `regulatory_genomics_toolbox` conda environment for footprinting analysis, motif analysis and differential motif analysis between the two conditions.
- `rnaseq`
    - `nf_rnaseq.sh`: Nextflow call used to preprocess the scRNA-seq data with the nfcore/rnaseq v1.4.2 pipeline.
    - `custom.config`: Custom configuration file for the Nextflow call specified in `nf_rnaseq.sh`.
    - `gene_level_summary.R`: Script to create a SummarizedExperiment object from the output files of the Nextflow pipeline.


