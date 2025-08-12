#!/bin/bash

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# indicate to run w/round 1 or round2
# sbatch  --export=round=1 ~/scripts/gather_applybqsr.slurm
module load apptainer/1.1.6
module load singularity/3.8.3
module load singularity/ce-3.11.0
module load gatk4/4.1.8.0--py38h37ae868_0

# define variable for bqsr round


wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

version='gatk --version'

echo "Run GATK GatherBQSRReports with ${version}  on  $(date) for round ${round}"

gatk GatherBQSRReports \
`(cat ${wd}/bqsr/bqsr${round}_reports.txt)` \
--output ${wd}/bqsr/recal${round}_data_combined.table

echo "Applying Recalibration table to bams on $(date) for round ${round}"

gatk ApplyBQSR \
-R ${ref} \
`(cat ${wd}/bqsr/bam_input.txt)` \
--bqsr-recal-file ${wd}/bqsr/recal${round}_data_combined.table \
-O ${wd}/bqsr/recal${round}_output.bam

mv applybqsr.* ${wd}/bqsr/std_outs/

