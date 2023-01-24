#!/bin/bash

#SBATCH -J dfpL                                 # Job name
#SBATCH -p cpu                              # Queue name
#SBATCH -N 1                                    # Number of nodes
#SBATCH -n 32                                   # Number of cores
#SBATCH --mem=128G                              # Memory pool for all cores
#SBATCH -t 3-00:00:00                          # Run time (hh:mm:ss)
#SBATCH -e dfpL.err                             # STDERR
#SBATCH -o dfpL.out                             # STDOUT
#SBATCH --mail-type=ALL                         # Type of notifications sent by e-mail
#SBATCH --mail-user=gregorio.alanis@crick.ac.uk       # e-mail address

# Load all necessary modules
ml Anaconda3

# Variable definitions
sample1=8cell.mRp.clN
prefix1=8cell

sample2=4cell.mRp.clN
prefix2=4cell

outDir=footprints
genome=hg38

# Activate the conda environment for RGT
source activate rgt_tools

# Remove peaks from non-standard contigs
awk '!/^chr.{1,2}_/' results/bwa/mergedReplicate/macs/narrowPeak/${sample1}_peaks.narrowPeak > ${outDir}/${sample1}.narrowPeak
awk '!/^chr.{1,2}_/' results/bwa/mergedReplicate/macs/narrowPeak/${sample2}_peaks.narrowPeak > ${outDir}/${sample2}.narrowPeak

# Identify footprints
rgt-hint footprinting --atac-seq --organism=$genome --output-location=$outDir --output-prefix=$prefix1 results/bwa/mergedReplicate/${sample1}.sorted.bam ${outDir}/${sample1}.narrowPeak 
rgt-hint footprinting --atac-seq --organism=$genome --output-location=$outDir --output-prefix=$prefix2 results/bwa/mergedReplicate/${sample2}.sorted.bam ${outDir}/${sample2}.narrowPeak

# Generate tracks for the genome browser
rgt-hint tracks --bc --bigWig --organism=$genome --output-location=$outDir results/bwa/mergedReplicate/${sample1}.sorted.bam ${outDir}/${sample1}.narrowPeak  --output-prefix=${prefix1}_bw
rgt-hint tracks --bc --bigWig --organism=$genome --output-location=$outDir results/bwa/mergedReplicate/${sample2}.sorted.bam ${outDir}/${sample2}.narrowPeak  --output-prefix=${prefix2}_bw

# Perform motif analysis
rgt-motifanalysis matching --organism=$genome --output-location=$outDir --filter "database:hocomoco,jaspar_vertebrates;species:sapiens" --input-files ${outDir}/${prefix1}.bed ${outDir}/${prefix2}.bed

# Perform differential motif analysis between the two conditions
rgt-hint differential --organism=$genome --bc --nc 30 --mpbs-files=${outDir}/${prefix1}_mpbs.bed,${outDir}/${prefix2}_mpbs.bed --reads-files=results/bwa/mergedReplicate/${sample1}.sorted.bam,results/bwa/mergedReplicate/${sample2}.sorted.bam --conditions=${prefix1},${prefix2} --output-location=$outDir

# Deactivate the conda environment
conda deactivate

