#!/bin/bash

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load apptainer/1.1.9
module load samtools/1.19.2
module load bcftools/1.19

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# Specify the path to the config file
config=${wd}/config_files/depth_filter_calcs.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
min_depth=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
max_depth=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the sample name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."

# use samtools to generate a consensus sequence

 version='bcftools --version'

echo "Generating consensus sequence for ${sample} with min DP ${min_depth} and max DP ${max_depth} filters on $(date) with ${version}"

mkdir -p ${wd}/psmc
mkdir -p ${wd}/psmc/std_outs

bcftools mpileup -C50 -q20 -Q20 -f ${ref} ${wd}/bam/${sample}_aln-pe_addRG.bam \
	-r  NC_066344.1,NC_066345.1,NC_066346.1,NC_066347.1,NC_066348.1,NC_066349.1,NC_066350.1,NC_066351.1,NC_066352.1,NC_066353.1,NC_066354.1,NC_066355.1,NC_066356.1,NC_066357.1,NC_066358.1,NC_066359
.1,NC_066360.1,NC_066361.1 \
	| bcftools call -c - \
	| vcfutils.pl vcf2fq -d ${min_depth} -D ${max_depth} | gzip > ${wd}/psmc/${sample}.diploid.fq.gz

echo "done generating consensus sequence on $(date)" 

# use psmc tools to generate input file from consensus sequence

# activate psmc environment
source activate ${conda_envs}/psmc

version='${software}/psmc/psmc --version'

echo "Generating psmc input file for ${sample} on $(date) with ${version}"

#${software}/psmc/utils/fq2psmcfa -q20  ${wd}/psmc/${sample}.diploid.fq.gz > ${wd}/psmc/${sample}.diploid.psmcfa

echo "done reformatting input file on $(date)" 

echo "running psmc on ${sample} on $(date)"

${software}/psmc/psmc -N30 -t15 -r5 -p "4+20*2+6" -o ${wd}/psmc/${sample}.psmc_output.psmc ${wd}/psmc/${sample}.diploid.psmcfa

echo "done running psmc on $(date)"
