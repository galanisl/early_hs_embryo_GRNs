#!/bin/sh
## LOAD REQUIRED MODULES
ml purge
ml Nextflow/19.10.0
ml Singularity/2.6.0-foss-2016b

## RUN PIPELINE
## SEE FULL PARAMETER LIST: https://github.com/nf-core/atacseq/blob/master/docs/usage.md 
nextflow run nf-core/atacseq \
    --genome GRCh38 \
    --input 'design_hs_atac.csv' \
    --single_end \
    --narrow_peak \
    --email YOUR_EMAIL_ADDRESS \
    -profile crick \
    -r 1.1.0 \
    -resume \
    -with-singularity /path/to/singularity/image/nfcore-atacseq-1.1.0.img

