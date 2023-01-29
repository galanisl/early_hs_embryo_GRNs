# Scripts for QC, filtering and normalisation of scRNA-seq and liATAC-seq data

The following scripts and files are provided:

- `00_scrnaseq_sample_selection_blastocyst.R`: Script to select high quality EPI, PE and TE samples from Blakeley *et al.*, Development 2015; Petropoulos *et al.*, Cell 2016 and Yan *et al.*, Nature Structural & Molecular Biology 2013.
- `01_scrnaseq_qc_processing_blastocyst.R`: Script for QC, filtering, integration and normalisation of EPI, PE and TE samples to generate the `SingleCellExperiment` object (i.e. the scRNA-seq data) used in this work.
- `02_scrnaseq_exprdata_to_matrix.R`: Script to export assays in the blastocyst `SingleCellExperiment` to matrix format.
- `01_licatseq_clean_CA.R`: Script to clean TF-gene associations from the liATAC-seq data and produce the TF-gene interactions used in this work.
- `gene_expr`:
    - `blakeley_salmon_fpkm_tpm.RDS`: Output of the `nf-core/rnaseq` processing pipeline for the Blakeley *et al.*, Development 2015 dataset.
    - `blakeley_sample_alias.tsv` Table with the different sample identifiers for the Blakeley *et al.*, Development 2015 dataset.
    - `blakeley_sample_details.csv`: Sample metadata for Blakeley *et al.*, Development 2015 dataset.
    - `petro_salmon_fpkm_tpm.RDS`: Output of the `nf-core/rnaseq` processing pipeline for the Petropoulos *et al.*, Cell 2016 dataset.
    - `petro_sample_details.csv`: Sample metadata for Petropoulos *et al.*, Cell 2016 dataset.
    - `yan_salmon_fpkm_tpm.RDS`: Output of the `nf-core/rnaseq` processing pipeline for the Yan *et al.*, Nature Structural & Molecular Biology 2013 dataset.
    - `yan_sample_alias.tsv`: Table with the different sample identifiers for the Yan *et al.*, Nature Structural & Molecular Biology 2013 dataset.
    - `yan_sample_details.csv`: Sample metadata for Yan *et al.*, Nature Structural & Molecular Biology 2013 dataset.
    - `samples_Wamaitha_NatComm2020.csv`: List of blastocyst samples used in Wamaitha *et al.*, Nature Communications 2020.
    - `stirparo_revlineages_lateBlast.csv`: Revised lineages for samples in Blakeley *et al.*, Development 2015; Petropoulos *et al.*, Cell 2016 and Yan *et al.*, Nature Structural & Molecular Biology 2013 according to Stirparo *et al.*, Development 2018.
- `licat_peaks`:
    - `ann_ICM_binding_sites.csv`: Regions of open chromatin in the ICM with TF motifs enriched and nearest promoters. This is the output of the `nf-core/atacseq` pipeline, followed by TF footprinting and motif enrichment analysis.
    - `ann_TE_binding_sites.csv`: Regions of open chromatin in the TE with TF motifs enriched and nearest promoters. This is the output of the `nf-core/atacseq` pipeline, followed by TF footprinting and motif enrichment analysis.
- `final_data`:
    - `sce_integrated_lateBlast.RData`: `SingleCellExperiment` experiment object with the scRNA-seq data before QC, filtering, integration, etc.
    - `sce_integrated_lateBlast_clean.RData`: `SingleCellExperiment` experiment object with the scRNA-seq data used in this work.
   - `icm_tf_binding.rds`: List of ICM-specific TF-gene interactions used in this work.
   - `te_tf_binding.rds`: List of TE-specific TF-gene interactions used in this work.
   - `regulators.rds`: List of TF with enriched motifs in the ICM and TE chromatin accessibility data. 