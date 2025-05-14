#!/bin/bash

#SBATCH -p atlas
#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="trimmomatic"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o trim.%a.%A.out                    # File to which stdout will be written
#SBATCH -e trim.%a.%A.err 
#SBATCH --account=beenome100

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda/4.12.0

source activate /project/beenome100/conda_envs/trimmomatic

wd=/90daydata/beenome100/rena_in_progress/${species}

# Specify the path to the config file
config=${wd}/config_files/trim_array.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the sample name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."

 version='trimmomatic PE --version'

echo "Generating trimmed output on $(date) with ${version}"

trimmomatic PE -threads 48 \
${wd}/raw_data//${sample}_R1.fastq.gz \
${wd}/raw_data/${sample}_R2.fastq.gz \
${wd}/trimmed_fastq/${sample}_R1.trimmed.fastq \
${wd}/trimmed_fastq/${sample}_R1un.trimmed.fastq \
${wd}/trimmed_fastq/${sample}_R2.trimmed.fastq \
${wd}/trimmed_fastq/${sample}_R2un.trimmed.fastq \
ILLUMINACLIP:adapters_iTru.fa:2:30:10:8:TRUE \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:5:20 \
MINLEN:36

mv trim.* ${wd}/trimmed_fastq/std_outs/

