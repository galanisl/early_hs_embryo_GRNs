#!/bin/sh
## LOAD REQUIRED MODULES
ml purge
ml Nextflow/19.10.0
ml Singularity/2.6.0-foss-2016b

## RUN PIPELINE
## SEE FULL PARAMETER LIST: https://github.com/nf-core/rnaseq/blob/master/docs/usage.md
nextflow run nf-core/rnaseq \
    --fasta /path/to//Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --gtf /path/to/Homo_sapiens.GRCh38.99.gtf \
    --star_index /path/to/star/index/if/available \
    --igenomes_base false \
    --singleEnd \
    --reads 'fastq/*.fastq.gz' \
    --pseudo_aligner salmon \
    --email YOUR_EMAIL_ADDRESS \
    -profile crick \
    -r 1.4.2 \
    -with-singularity /path/to/singularity/image/nfcore-rnaseq-1.4.2.img \
    -c custom.config \
    --skipMultiQC \
    -resume 

