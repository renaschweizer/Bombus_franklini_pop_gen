#!/bin/bash


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load apptainer/1.1.6
module load singularity/3.8.3
module load singularity/ce-3.11.0
module load gatk4/4.1.8.0--py38h37ae868_0

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# Specify the path to the config file
config=${wd}/config_files/genome_contigs_array.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
contig=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

version=`gatk --version`

echo "Run GATK BaseRecalibrator ROUND 2 ${version} for ${contig} on  $(date)"

gatk BaseRecalibrator \
--input ${wd}/bqsr/recal1_output.bam \
--reference ${ref} \
--known-sites ${wd}/mpileup_gt/${species}_known_sites.vcf \
--intervals ${contig} \
--output ${wd}/bqsr/recal2_data_${contig}.table

mv bqsr2.* ${wd}/bqsr/std_outs/

