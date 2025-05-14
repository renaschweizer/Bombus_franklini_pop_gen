#!/bin/bash

#SBATCH -p atlas
#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="addRG"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o addRG.%a.%A.out                    # File to which stdout will be written
#SBATCH -e addRG.%a.%A.err 
#SBATCH --account=beenome100

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load apptainer/1.1.6
module load singularity/3.8.3
module load singularity/ce-3.11.0
module load picard/2.23.2--0

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# Specify the path to the config file
config=${wd}/config_files/addRG_info_array.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
lane_num=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
run_id=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the sample name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."
echo "lane num: ${lane_num} and run id: ${run_id}"

 version='picard --version'

echo "Adding read group information to bam files on $(date) with ${version}"

## DEFINE READ GROUP INFO 

## RGID=String	Read Group ID Default value: 1. used by GATK
RGID_specify=${run_id}_${lane_num}
## RGLB=String	Read Group Library Required. used by GATK/picardmarkdups
RGLB_specify=${sample} 
## RGPL=String	Read Group platform (e.g. illumina, solid) Required. used by GATK
RGPL_specify=ILLUMINA
## RGPU=String	Read Group platform unit (eg. run barcode of run Id + lane number) Required. 
RGPU_specify=${run_id}_${lane_num}
## RGSM=String	Read Group sample name Required. used by GATK
RGSM_specify=${sample}
####################
#	addReadGroupBam
#	adds read group information into the bam file
# 	check with samtools view -H /path/to/my.bam | grep '^@RG'
####################

mkdir ${wd}/bam/temp
 
picard AddOrReplaceReadGroups \
INPUT=${wd}/bam/${sample}_aln-pe.s.u.bam \
OUTPUT=${wd}/bam/${sample}_aln-pe_addRG.bam \
SORT_ORDER=coordinate \
RGID=${RGID_specify} \
RGLB=${RGLB_specify} \
RGPL=${RGPL_specify} \
RGPU=${RGPU_specify} \
RGSM=${RGSM_specify} \
CREATE_INDEX=True \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=${wd}/bam/temp

mv addRG.* ${wd}/bam/std_outs/

