#!/bin/bash

#SBATCH --time=04:00:00   
#SBATCH --nodes=1   
#SBATCH --ntasks-per-node=48   
#SBATCH --job-name="FastQC_clean"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o FastQC.%a.%A.out                    
#SBATCH -e FastQC.%a.%A.err 


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load fastqc/0.11.9

wd=/90daydata/beenome100/rena_in_progress/${species}

# Specify the path to the config file to run fastqc on raw or trimmed reads
#config=${wd}/config_files/raw_fastq_array.txt
config=${wd}/config_files/fastq_trimmed_array.txt


# Extract the file name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, the same name, and the sex of the sample
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."

 version='fastqc --version'

echo "Generating fastqc output  on $(date) with ${version}"

fastqc -o ${wd}/fastqc_reports -t 2 --noextract ${wd}/${sample} 

mv FastQC.* ${wd}/fastqc_reports/std_outs/

